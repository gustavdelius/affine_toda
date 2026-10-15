# Phase-2 referee report: Part III (chapters 9–11)

- **Date:** 2026-10-15
- **Model:** Gemini 4 Argon (High)
- **Commit reviewed:** `a03edd4896c927253721352a5bb253521c7f0f39`

---

## Item 1. The $\Theta$ selection rule and its entry in Appendix E (`11-reality.qmd`, `appendices/E-original-results.qmd:8`)

### Overlap

- **Overlap with Phase 1:** **in both (partial)**.
  - **Covered in Phase 1:** Rank 1 (`11-reality.qmd:32–41, 47`, `10-solitons.qmd:180`) and Rank 2 (`11-reality.qmd:43–49`, `appendices/E-original-results.qmd:8`). Phase 1 checked the action of $\Theta:\vec\phi\mapsto -\bar{\vec\phi}$ (`\tau_j\mapsto\bar\tau_j`) on Hirota data across all simply-laced algebras, identified the hidden hypotheses (additivity/independence from $\xi_i$, $\operatorname{Re}\cosh\theta_i>0$, $\prod_{k<l}A_{a_k a_l}\notin\{0,\infty\}$, and a $\Theta$-eigenbasis when an exponent is degenerate), noted the missing `::: {#prp-...}` statement div and status label in `11-reality.qmd` versus `appendices/E-original-results.qmd:8`, checked Khastgir–Sasaki (`hep-th/9507001`) on complex-energy vs negative-energy solutions, and symbolically constructed and numerically integrated the spin-2 conserved current of $a_2^{(1)}$ on single solitons and breathers.
  - **New in Phase 2:**
    1. Checking the spin-4 conserved current of $a_2^{(1)}$ and Freeman's general soliton–antisoliton parity formula (`hep-th/9408092`, eqs. (3.48)–(3.49)) to determine $\epsilon_Q$ for all spins $s$.
    2. Prior-art comparison against the four papers of Fring and collaborators (`arXiv:1602.05465`, `arXiv:1608.01691`, `arXiv:2007.15425`, `arXiv:2110.06825`) on $\mathcal{PT}$-symmetric complex solitons with real energies.
    3. Direct comparison with gap H6 in `part9-foundations/44-twelve-gaps.qmd:33, 81`.

### Check

1. **Action of $\Theta$ on Hirota data (`eq-theta-data`, `11-reality.qmd:20–24`) and hidden hypotheses (`11-reality.qmd:43–45`):**
   - For any affine Toda theory, conjugating the imaginary-coupling equation of motion (`eq-toda-eom`, $\partial_+\partial_-\vec\phi = \frac{im^2}{\beta}\sum_j n_j\vec\alpha_j e^{i\beta\vec\alpha_j\cdot\vec\phi}$) gives $\partial_+\partial_-(-\bar{\vec\phi}) = \frac{im^2}{\beta}\sum_j n_j\vec\alpha_j e^{i\beta\vec\alpha_j\cdot(-\bar{\vec\phi})}$, so $\Theta:\vec\phi\mapsto -\bar{\vec\phi}$ maps solutions to solutions.
   - On the Hirota ansatz $\vec\phi = \frac{i}{\beta}\sum_j \frac{2\vec\alpha_j}{\alpha_j^2}\ln\tau_j$ (`eq-hirota-ansatz`), since $\beta>0$ and $\vec\alpha_j\in\mathbb R^r$, we have $-\bar{\vec\phi} = \frac{i}{\beta}\sum_j \frac{2\vec\alpha_j}{\alpha_j^2}\ln\bar\tau_j$ (up to branch cuts of $\ln\tau_j$, which shift $\vec\phi$ by a vacuum in $\frac{2\pi}{\beta}\Lambda_W^\vee$). Thus $\Theta$ acts by $\tau_j \mapsto \bar\tau_j$.
   - For $a_n^{(1)}$ (`eq-n-soliton`), $\overline{\omega^{ja}} = \omega^{-ja} = \omega^{j(h-a)}$, $\overline{E_k(\theta_k,\xi_k)} = E_k(\bar\theta_k,\bar\xi_k)$, and from `eq-interaction-coefficient`:
     $$\overline{A_{ab}(\theta)} = \frac{\cosh\bar\theta - \cos\frac{\pi(a-b)}{h}}{\cosh\bar\theta - \cos\frac{\pi(a+b)}{h}} = \frac{\cosh\bar\theta - \cos\frac{\pi((h-a)-(h-b))}{h}}{\cosh\bar\theta - \cos\left(2\pi - \frac{\pi((h-a)+(h-b))}{h}\right)} = A_{h-a,h-b}(\bar\theta).$$
     Thus $\bar\tau_j$ is the $N$-soliton solution with constituent data $\{(\bar a_k,\bar\theta_k,\bar\xi_k)\}$ where $\bar a_k = h - a_k$. (As verified in Phase 1 Rank 2 using MacKay–McGhee `hep-th/9208057`, $\tau_j\mapsto\bar\tau_j$ likewise maps $(a,\theta,\xi)\mapsto(\bar a,\bar\theta,\bar\xi)$ in $d_n^{(1)}$ and $e_{6,7,8}^{(1)}$.)
   - **Hidden assumptions in the selection rule (`11-reality.qmd:45`):**
     If $Q = \int_{-\infty}^\infty J_Q(\vec\phi)\,dx$ is a local conserved charge whose density satisfies $J_Q(-\bar{\vec\phi}) = \epsilon_Q\overline{J_Q(\vec\phi)}$ ($\epsilon_Q = \pm 1$), then
     $$Q[\Theta\vec\phi] = \epsilon_Q\overline{Q[\vec\phi]}.$$
     However, $\Theta$-invariance of the asymptotic multiset $\{(a_i,\theta_i)\}$ does **not** imply $\Theta\vec\phi = \vec\phi$ as a field, because $\Theta$ maps $\xi_i \mapsto \bar\xi_i$, which need not equal $\xi_{\pi(i)}$ (indeed, for a regular $a_n^{(1)}$ breather, $\xi_2 = \bar\xi_1 + i\pi$, so $\bar\xi_1 = \xi_2 - i\pi \ne \xi_2$; and under a real spacetime translation $(x,t)\mapsto(x-x_0,t-t_0)$ on a pair with $\theta_2 = \bar\theta_1$, $\delta\xi_2 = \overline{\delta\xi_1}$, so $\xi_2 - \bar\xi_1$ is strictly invariant under real translations and cannot be removed by shifting the origin). Therefore, concluding $Q[\Theta\vec\phi] = Q[\vec\phi]$ requires the assumption that **$Q$ depends only on the constituent multiset $\{(a_i,\theta_i)\}$ and is independent of the position/phase parameters $\xi_i$** (which follows from additivity $Q = \sum_i Q_1(a_i,\theta_i)$ when $\operatorname{Re}\cosh\theta_i > 0$, $\prod_{k<l}A_{a_k a_l}(\theta_k-\theta_l)\notin\{0,\infty\}$, and no zero of $\tau_j$ lies on the integration contour). At `11-reality.qmd:45`, the first bullet states $Q = \epsilon_Q\bar Q$ unconditionally, mentioning additivity only in the second bullet (`11-reality.qmd:46`) and in the status label of `appendices/E-original-results.qmd:8`.

2. **Check of $\epsilon_Q = \pm 1$ for spin 1, spin 2, and spin 4 in $a_2^{(1)}$ (`referee_scratch/classical_09_11/check_ch11_reality_breathers.py` Check 4 and `check_phase2_items.py` Item 1):**
   - In terms of $\vec\psi = i\beta\vec\phi = -\sum_j \vec\alpha_j^\vee\ln\tau_j$, the imaginary-coupling field equation (`eq-toda-eom`) is $\partial_+\partial_-\vec\psi = -m^2\sum_j n_j\vec\alpha_j e^{\vec\alpha_j\cdot\vec\psi}$, which has strictly real coefficients, and $\Theta:\vec\phi\mapsto -\bar{\vec\phi}$ is plain complex conjugation $\vec\psi\mapsto\bar{\vec\psi}$.
   - On a static or real-rapidity single soliton of species $a$ at $\xi=0$, $\overline{\tau_j(a)} = \tau_j(\bar a)$, so $\Theta$ coincides with charge conjugation $C: a\leftrightarrow\bar a$ (the Dynkin diagram involution permuting conjugate nodes $j\leftrightarrow\bar j$).
   - By Freeman (`hep-th/9408092`, eqs. (3.48)–(3.49), quoted in Item 3 below), for any untwisted affine Toda theory at spin $N$, the single-soliton charges form a left eigenvector of the Cartan matrix of $g$ with eigenvalue $2(1-\cos(N\pi/h)) = 4\sin^2\frac{N\pi}{2h}$, and obey the soliton–antisoliton relation
     $$Q_{-N}(\bar a) = (-1)^{N-1} Q_{-N}(a).$$
     When $J_{Q_N}$ is normalized so that the single-soliton charges $q_a^{(N)}$ are real (as in Freeman line 704 and `11-reality.qmd:46`), evaluating $Q_N[\Theta\vec\phi_a] = \epsilon_{Q_N}\overline{Q_N[\vec\phi_a]}$ on a real-rapidity 1-soliton gives $q_{\bar a}^{(N)} = \epsilon_{Q_N} q_a^{(N)}$, so
     $$\epsilon_{Q_N} = \frac{q_{\bar a}^{(N)}}{q_a^{(N)}} = (-1)^{N-1}$$
     for every spin $N$ (in a $\Theta$/$C$-eigenbasis when an exponent has multiplicity $>1$, as in $d_{2m}^{(1)}$ at $N=2m-1$).
   - Specifically for $a_2^{(1)}$ ($h=3$, exponents $N = 1, 2, 4, 5, \dots \not\equiv 0\pmod 3$), $q_a^{(N)} \propto \sin\frac{N\pi a}{3}$:
     - **Spin $N=1$ ($E, P$):** $T_2 = \frac12(\partial_+\vec\phi)^2$ has real coefficients and even degree $2$ in $\vec\phi$, so $J_{Q_1}(-\bar{\vec\phi}) = +\overline{J_{Q_1}(\vec\phi)}$ ($\epsilon_{Q_1} = +1$), and $q_2^{(1)}/q_1^{(1)} = \sin(2\pi/3)/\sin(\pi/3) = +1$ (Cartan eigenvalue $1$).
     - **Spin $N=2$:** In the $\mathbb Z_2$-conjugation basis $\phi_1 = \frac{1}{\sqrt{2}}(-\vec\alpha_0)\cdot\vec\phi$ ($C$-even), $\phi_2 = \frac{1}{\sqrt{6}}(\vec\alpha_1-\vec\alpha_2)\cdot\vec\phi$ ($C$-odd), the conserved current $\partial_- J_{+3} + \partial_+ J_{-1} = 0$ is (`check_ch11_reality_breathers.py`):
       $$J_{+3} = \partial_+\phi_2\big((\partial_+\phi_2)^2 - 3(\partial_+\phi_1)^2\big) - \frac{3i}{\sqrt{2}\beta}\big(\partial_+\phi_1\,\partial_+^2\phi_2 - \partial_+^2\phi_1\,\partial_+\phi_2\big),$$
       $$J_{-1} = \frac{3m^2}{2\beta^2}\Big[2e^{i\beta\vec\alpha_0\cdot\vec\phi}\partial_+\phi_2 + e^{i\beta\vec\alpha_1\cdot\vec\phi}(\sqrt{3}\partial_+\phi_1 - \partial_+\phi_2) - e^{i\beta\vec\alpha_2\cdot\vec\phi}(\sqrt{3}\partial_+\phi_1 + \partial_+\phi_2)\Big].$$
       Under $\Theta:\vec\phi\mapsto -\bar{\vec\phi}$, the cubic term in $J_{+3}$ has real coefficients and odd degree $3$, the quadratic derivative term has a purely imaginary coefficient $-3i/(\sqrt{2}\beta)$ and even degree $2$, and $J_{-1}$ has real coefficients times $e^{i\beta\vec\alpha_j\cdot\vec\phi}\partial_+\phi_k$. Hence $J_{Q_2}(-\bar{\vec\phi}) = -\overline{J_{Q_2}(\vec\phi)}$ ($\epsilon_{Q_2} = -1$), and $q_2^{(2)}/q_1^{(2)} = \sin(4\pi/3)/\sin(2\pi/3) = -1$ (Cartan eigenvalue $3$).
     - **Spin $N=4$:** Although $N=4 \equiv 1\pmod 3$ comes from the same modulo-$h$ exponent $\nu=1$ as $N=1$, in Freeman's eq. (3.48) the factor $\exp(\pm(N-\nu)\pi i/(2h)) = \exp(\pm i\pi/2) = \pm i$ flips the relative sign between white node $1$ and black node $2$ and introduces a factor of $i$. In the $W_3$ algebra, the spin-4 density at weight $5$ is proportional to the composite current $:T_2 W_3:$ (product of $C$-even spin-2 $T_2$ with $\epsilon=+1$ and $C$-odd spin-3 $W_3$ with $\epsilon=-1$), so with real single-soliton normalization it satisfies $J_{Q_4}(-\bar{\vec\phi}) = -\overline{J_{Q_4}(\vec\phi)}$ ($\epsilon_{Q_4} = (-1)^{4-1} = -1$), with $q_2^{(4)}/q_1^{(4)} = \sin(8\pi/3)/\sin(4\pi/3) = -1$ (Cartan eigenvalue $3$). Output of `check_phase2_items.py`:
       ```text
       a_2^(1) spin N=1: Cartan eigenvalue = 1.0, q_2/q_1 = +1.0, eps_{Q_1} = (-1)^(1-1) = +1
       a_2^(1) spin N=2: Cartan eigenvalue = 3.0, q_2/q_1 = -1.0, eps_{Q_2} = (-1)^(2-1) = -1
       a_2^(1) spin N=4: Cartan eigenvalue = 3.0, q_2/q_1 = -1.0, eps_{Q_4} = (-1)^(4-1) = -1
       a_2^(1) spin N=5: Cartan eigenvalue = 1.0, q_2/q_1 = +1.0, eps_{Q_5} = (-1)^(5-1) = +1
       ```

3. **Comparison with prior art by Fring and collaborators (`11-reality.qmd:96`, `appendices/E-original-results.qmd:8`):**
   We read all four papers:
   - **Cen and Fring (`arXiv:1602.05465`, *J. Phys. A* 49 (2016) 365202, "Complex solitons with real energies"):**
     In Section 2.1.3 and Section 2.3.1 (eqs. (50)–(60), complex sine-Gordon), Cen and Fring consider complex-shifted $N$-soliton solutions and state two mechanisms for reality of the energy and higher integrals of motion:
     1. *Single solitons:* the field equation is invariant under the antilinear $\mathcal{PT}$ symmetry $\mathcal{PT}_-: x\mapsto -x, t\mapsto -t, i\mapsto -i, \phi\mapsto -\phi$ (whose field part $\phi\mapsto -\bar\phi$ is $\Theta$). Even when a complex shift $\mu$ breaks $\mathcal{PT}$ of $\phi(x,t)$, shifting the real part of the soliton centre restores $\mathcal{PT}$ up to a topological vacuum shift, which makes the energy density satisfies $\mathcal H(x,t) = \mathcal H^*(-x,-t)$ and forces the integrated 1-soliton energy to be real.
     2. *Multi-solitons:* when $N \ge 2$ solitons have distinct real rapidities and different imaginary shifts $\mu_j$, the $N$-soliton solution is **not** $\mathcal{PT}$-symmetric as a field at finite time (`arXiv:1602.05465`, Section 2.1.3: *"We cannot restore the $\mathcal{PT}$-symmetry by a simple shift in $x$ cleanly... Thus we can not argue directly that the energy should be real"*). Instead, Cen and Fring argue that because the $N$-soliton solution separates asymptotically as $t\to\pm\infty$ into $N$ single solitons with the same conserved total energy, the $N$-soliton energy is the sum of the real single-soliton energies (`arXiv:1602.05465`, Section 2.1.3 and eq. (60)).
   - **Cen, Correa and Fring (`arXiv:1608.01691`, *J. Math. Phys.* 58 (2017) 032901, "Time-delay and reality conditions for complex solitons"):**
     In Section 3 ("Reality of All Charges from Asymptotic $\mathcal{PT}$-symmetry", eqs. (70)–(78)), Cen, Correa and Fring generalize this to **all higher conserved charges $I_n$** of $N$-soliton solutions:
     > *"Here we generalize this construction iteratively to guarantee the reality of all charges expressions $I_n$. A central property of all gradient density functions $w_n$ related to conserved quantities directly is that they transforms with a definite parity under a $\mathcal{PT}$-symmetry transformation $\mathcal{PT}(w_n)=(-1)^n w_n$."* (`arXiv:1608.01691`, p. 11)
     > *"Thus by standard arguments for integrable systems, we conclude that the total charge of any $N$-soliton solution may be decomposed into a sum of 1-soliton charges $I_n(u^{(N)}) = \sum_{k=1}^N I_n(u_k^{(1)})$."* (`arXiv:1608.01691`, eq. (77))
   - **Fring and Taira (`arXiv:2007.15425`, *J. Phys. A* 53 (2020) 455701, "Complex BPS solitons with real energies from duality"):**
     Treats two-component complex scalar models with $\mathcal{PT}$ symmetry where topological charges (Bogomolny bounds) of single complex solitons are real when the vacua invaded at $x\to\pm\infty$ are $\mathcal{PT}$-invariant or rotated by a dual symmetry.
   - **Correa, Fring and Taira (`arXiv:2110.06825`, *Front. Phys.* 10 (2022) 821861, "Linearly stable and unstable complex soliton solutions with real energies in the Bullough-Dodd model"):**
     Studies imaginary-coupling Bullough–Dodd ($a_2^{(2)}$) with $\varphi = -i\beta\phi = \ln(\tau_0/\tau_1)$. Under $\mathcal{PT}_1: x\mapsto -x, t\mapsto -t, i\mapsto -i, \varphi\mapsto\varphi$ (which in terms of $\phi = i\varphi/\beta$ is $\phi\mapsto -\bar\phi$, i.e. $\Theta$ combined with spacetime reflection), they construct complex 1-solitons, 2-solitons, and breathers via Hirota's method and compute their real energies from the asymptotic total-derivative formula $E = [\ldots]_{-\infty}^\infty$.
   - **What is in Fring et al. vs what is genuinely new in Chapter 11:**
     - **Already in Fring et al. (`arXiv:1602.05465`, `arXiv:1608.01691`, `arXiv:2110.06825`):** (i) Using the antilinear symmetry $\vec\phi\mapsto -\bar{\vec\phi}$ (combined with $(x,t)\mapsto(-x,-t)$ as $\mathcal{PT}$) together with asymptotic additivity over constituent solitons to explain why multi-soliton configurations that are *not* $\mathcal{PT}$-symmetric as fields at finite time nonetheless have real energy and higher charges; and (ii) applying this to sine-Gordon (`arXiv:1602.05465`) and Bullough–Dodd $a_2^{(2)}$ (`arXiv:2110.06825`).
     - **Genuinely new in Chapter 11:** Cen and Fring (`arXiv:1602.05465`, Section 3) explicitly restricted their general reality mechanism to **real dispersion parameters/rapidities** (*"Here we have not allowed complex dispersion relations, i.e. keeping our parameters $\alpha,\beta$ real"*) and to single-field models where $\bar a = a$. What is new in Chapter 11 is:
       1. Formulating the selection rule directly on the **multiset of constituent Hirota data $\{(a_i,\theta_i)\}$** under $\Theta:(a,\theta,\xi)\mapsto(\bar a,\bar\theta,\bar\xi)$ (`eq-theta-data`) for general multicomponent affine Toda theories $\hat g$ (where conjugate species $\bar a \ne a$ exist);
       2. Extending the selection rule to **complex-rapidity conjugate pairs** $\{(a,\theta_0+iu),(\bar a,\theta_0-iu)\}$ (breathers) and using it to exclude unpaired complex-rapidity solutions (the Khastgir–Sasaki complex-energy solitons); and
       3. Distinguishing the weaker condition $\{(M_{a_i},\theta_i)\} = \{(M_{a_i},\bar\theta_i)\}$ for $C$-even charges ($\epsilon_Q = +1$, such as $E$ and $P$ on a single non-self-conjugate soliton $a\ne\bar a$ at real $\theta$) from full $\Theta$-invariance $\{(a_i,\theta_i)\} = \{(\bar a_i,\bar\theta_i)\}$ required for $C$-odd charges ($\epsilon_Q = -1$, such as the spin-2 and spin-4 charges in $a_2^{(1)}$, which are purely imaginary on a breather and real non-zero on a single soliton).
     - Consequently, `11-reality.qmd:96` is **wrong** to omit Cen–Fring (`arXiv:1602.05465`), Cen–Correa–Fring (`arXiv:1608.01691`), and Correa–Fring–Taira (`arXiv:2110.06825`) and to claim unqualifiedly that *"The formulation as a selection rule from $\Theta$ is new here"*. Both `11-reality.qmd:96` and `appendices/E-original-results.qmd:8` must cite those papers and state precisely what is new here.

4. **Check against gap H6 (`part9-foundations/44-twelve-gaps.qmd:33, 81` and `11-reality.qmd:49`):**
   - `44-twelve-gaps.qmd:33` defines **H6 (Energy reality and solution selection)** as:
     > *"Reality of soliton energies has been shown for specific solution classes. There are also complex-energy and negative-energy solutions [@khastgir1996], and singular ones. No general criterion distinguishes the solutions that should count as physical."*
   - `44-twelve-gaps.qmd:81` states in the summary table: `H6 | Classical selection rule ($\Theta$-invariance) | @sec-classical | closed: @prp-hirota-energy and the $\Theta$ selection rule`.
   - `11-reality.qmd:49` states: *"This closes gap H6 of @sec-critique, apart from the question of singularities."*
   - **Discrepancy:** As verified in Phase 1 Rank 1 by reading Khastgir and Sasaki (`hep-th/9507001`), Khastgir and Sasaki exhibit **two** distinct classes of unphysical solutions:
     1. Single-soliton Hirota solutions with complex rapidity $\theta = \theta_R + i\theta_I$ (`hep-th/9507001`, Section 3.2, eq. (3.12)), which have complex energy $M_a\cosh(\theta_R + i\theta_I)$ (with $\operatorname{Re}E > 0$ when $\cos\theta_I > 0$, or real negative energy when $\theta_I = \pi$ where $\operatorname{Re}\cosh\theta < 0$ and $\tau_j$ has real zeros). These are excluded by the $\Theta$ selection rule on Hirota data (together with $\operatorname{Re}\cosh\theta_i > 0$).
     2. **Negative-energy Cauchy solutions** (`hep-th/9507001`, Section 5.2, eqs. (5.12)–(5.13)): initial data $\vec\phi(x,0) = 0$ and purely imaginary initial velocity $\partial_t\vec\phi(x,0) = i\vec f(x)$ with real $\vec f\in L^2(\mathbb R,\mathbb R^r)$, for which $V_{\rm imag}(\vec\phi(x,0)) = 0$ and
        $$E = \int_{-\infty}^\infty \tfrac12\big(\partial_t\vec\phi(x,0)\big)^2\,dx = -\frac12\int_{-\infty}^\infty \vec f(x)^2\,dx < 0.$$
        Crucially, this solution has $\vec\phi(x,t)$ purely imaginary for all $(x,t)$, so it is **strictly $\Theta$-invariant as a field** ($\Theta\vec\phi = -\bar{\vec\phi} = \vec\phi$)! Thus $\Theta$-invariance alone does **not** exclude Khastgir and Sasaki's negative-energy radiation solutions; what excludes them is restricting to the Hirota multi-soliton sector with $\operatorname{Re}\cosh\theta_i > 0$ (`@prp-hirota-energy`), where $M_a > 0$ guarantees $E > 0$. This must be stated explicitly in `11-reality.qmd:47–49` and `44-twelve-gaps.qmd:33, 81`.

### Verdict

- **Verdict:** **wrong** (requires: (1) adding a `::: {#prp-theta-selection-rule}` statement div with status label `**[Theorem, given additivity of the higher charges]**` and proof in `11-reality.qmd:43–46` so that `appendices/E-original-results.qmd:8` cites a statement id as required by `AGENTS.md` and `check-refs.py`; (2) stating the hidden hypotheses of additivity/independence from $\xi_i$ and $\operatorname{Re}\cosh\theta_i > 0$ inside the statement; (3) citing Cen–Fring `arXiv:1602.05465`, Cen–Correa–Fring `arXiv:1608.01691`, and Correa–Fring–Taira `arXiv:2110.06825` in `11-reality.qmd:96` and narrowing the originality claim in `11-reality.qmd:96` and `appendices/E-original-results.qmd:8`; and (4) accurately distinguishing Khastgir–Sasaki's complex-energy soliton solutions from their $\Theta$-invariant negative-energy Cauchy solutions in `11-reality.qmd:47` and `10-solitons.qmd:180`).

### Paste-ready replacement text

- In `part3-classical/11-reality.qmd`, lines 43–50, replace:
  ```markdown
  On Hirota data, $\Theta$ acts as $\tau_j\mapsto\bar\tau_j$. That is, $(a,\theta,\xi)\mapsto(\bar a,\bar\theta,\bar\xi)$, where $\bar a$ is the conjugate species ($\omega^{ja}\mapsto\omega^{-ja}$). Consequently:

  - **Selection rule.** If the asymptotic data are $\Theta$-invariant, i.e. the multiset $\{(a_i,\theta_i)\}$ is mapped to itself by $(a,\theta)\mapsto(\bar a,\bar\theta)$, then every conserved charge satisfies $Q=\epsilon_Q\bar Q$, where $J_Q(-\bar{\vec\phi})=\epsilon_Q\overline{J_Q(\vec\phi)}$ defines $\epsilon_Q=\pm1$; in particular $E$ and $P$ are real. Breathers are of this kind.
  - **A weaker condition for $E$ and $P$.** By @prp-hirota-energy, $E$ and $P$ depend on the data only through the multiset $\{(M_{a_i},\theta_i)\}$. Since $M_{\bar a}=M_a$, they are real whenever this multiset is mapped to itself by $\theta\mapsto\bar\theta$, whether or not the species are conjugated as well. This covers every solution with real rapidities, for instance a single soliton with $\bar a\ne a$, which is not $\Theta$-invariant. Charges that are odd under charge conjugation $a\leftrightarrow\bar a$ need full $\Theta$-invariance. In $a_2^{(1)}$ the spin-2 charge of a single soliton is $q_ae^{2\theta}$, where by Freeman's computation [@freeman1995] the $q_a$ form the eigenvector of the $A_2$ Cartan matrix with eigenvalue $3$, so $q_2=-q_1$, and they can be normalized to be real. If the charge is additive over the constituents, as $E$ and $P$ are, a breather $\{(1,\theta_0+iu),(2,\theta_0-iu)\}$ has spin-2 charge $q_1\big(e^{2(\theta_0+iu)}-e^{2(\theta_0-iu)}\big)=2iq_1e^{2\theta_0}\sin2u$, which is imaginary. This agrees with the selection rule with $\epsilon_Q=-1$. Freeman treats single solitons only, so additivity at complex rapidities is assumed here.
  - **Complex energies.** The complex-energy solutions of Khastgir and Sasaki have asymptotic data that are not $\Theta$-invariant. They are not in the physical sector, for the same reason that complex eigenvalues of a $\Theta$-symmetric operator are not $\Theta$-invariant states (@sec-antilinear).

  This closes gap H6 of @sec-critique, apart from the question of singularities.
  ```
  with:
  ```markdown
  On Hirota data, $\Theta$ acts as $\tau_j\mapsto\bar\tau_j$. That is, $(a,\theta,\xi)\mapsto(\bar a,\bar\theta,\bar\xi)$, where $\bar a$ is the conjugate species ($\omega^{ja}\mapsto\omega^{-ja}$).

  ::: {#prp-theta-selection-rule}
  **[Theorem, given additivity of the higher charges]** Let $Q=\int dx\,J_Q(\vec\phi)$ be a local conserved charge whose density has definite $\Theta$-parity, $J_Q(-\bar{\vec\phi})=\epsilon_Q\overline{J_Q(\vec\phi)}$ with $\epsilon_Q=\pm1$, and assume that on regular Hirota multi-soliton solutions (with $\operatorname{Re}\cosh\theta_i>0$ and $\prod_{k<l}A_{a_ka_l}(\theta_k-\theta_l)\notin\{0,\infty\}$) $Q$ is additive over the constituents, $Q=\sum_i Q_1(a_i,\theta_i)$, and in particular independent of the position parameters $\xi_i$. If the asymptotic data $\{(a_i,\theta_i)\}$ are $\Theta$-invariant, i.e. mapped to themselves as a multiset by $(a,\theta)\mapsto(\bar a,\bar\theta)$, then $Q=\epsilon_Q\bar Q$; in particular $E$ and $P$ ($\epsilon_E=\epsilon_P=+1$) are real and positive. Breathers are of this kind.
  :::

  ::: {.proof}
  Conjugating the $\tau$-functions maps a Hirota multi-soliton solution $\vec\phi$ with parameters $\{(a_i,\theta_i,\xi_i)\}$ to the solution $\Theta\vec\phi=-\bar{\vec\phi}$ with parameters $\{(\bar a_i,\bar\theta_i,\bar\xi_i)\}$. Integrating $J_Q(\Theta\vec\phi)=\epsilon_Q\overline{J_Q(\vec\phi)}$ gives $Q[\Theta\vec\phi]=\epsilon_Q\overline{Q[\vec\phi]}$. By additivity and independence of the $\xi_i$, $\Theta$-invariance of the multiset $\{(a_i,\theta_i)\}$ implies $Q[\Theta\vec\phi]=Q[\vec\phi]$, so $Q=\epsilon_Q\bar Q$. For $E$ and $P$, additivity and positivity of $E$ when $\operatorname{Re}\cosh\theta_i>0$ follow from @prp-hirota-energy.
  :::

  Two consequences should be noted:

  - **A weaker condition for $E$ and $P$.** By @prp-hirota-energy, $E$ and $P$ depend on the data only through the multiset $\{(M_{a_i},\theta_i)\}$. Since $M_{\bar a}=M_a$, they are real whenever this multiset is mapped to itself by $\theta\mapsto\bar\theta$, whether or not the species are conjugated as well. This covers every solution with real rapidities, for instance a single soliton with $\bar a\ne a$, which is not $\Theta$-invariant. Charges that are odd under charge conjugation $a\leftrightarrow\bar a$ need full $\Theta$-invariance. By Freeman's computation [@freeman1995, eq. 3.49], at spin $s$ the single-soliton charges obey $q_{\bar a}=(-1)^{s-1}q_a$ and can be normalized to be real, which fixes $\epsilon_Q=(-1)^{s-1}$. In $a_2^{(1)}$ the spin-2 charge of a single soliton is $q_ae^{2\theta}$, where the $q_a$ form the eigenvector of the $A_2$ Cartan matrix with eigenvalue $3$, so $q_2=-q_1$ (and similarly at spin $4$). If the charge is additive over the constituents, as $E$ and $P$ are, a breather $\{(1,\theta_0+iu),(2,\theta_0-iu)\}$ has spin-2 charge $q_1\big(e^{2(\theta_0+iu)}-e^{2(\theta_0-iu)}\big)=2iq_1e^{2\theta_0}\sin2u$, which is imaginary. This agrees with @prp-theta-selection-rule with $\epsilon_Q=-1$. Freeman treats single solitons only, so additivity at complex rapidities is assumed here.
  - **Complex and negative energies.** The complex-energy single-soliton solutions of Khastgir and Sasaki [@khastgir1996] have complex rapidities whose asymptotic data are not $\Theta$-invariant; they are excluded for the same reason that complex eigenvalues of a $\Theta$-symmetric operator do not belong to $\Theta$-invariant states (@sec-antilinear). (Their negative-energy examples with $\operatorname{Re}\cosh\theta<0$ violate the decay condition $\operatorname{Re}\cosh\theta_i>0$ and are singular on the real line, while their smooth negative-energy Cauchy data with $\vec\phi(x,0)=0$ and imaginary $\partial_t\vec\phi(x,0)$ are $\Theta$-invariant radiation configurations outside the Hirota soliton sector.)

  Together with @prp-hirota-energy, this closes gap H6 of @sec-critique for Hirota soliton solutions, apart from the question of singularities.
  ```
- In `part3-classical/11-reality.qmd`, line 96 (`#sec-classical-notes`), replace:
  ```markdown
  Reality of soliton energies was shown by Olive, Turok and Underwood [-@olive1993a; -@olive1993b] and Zhu and Caldi [@zhu1995], and of the higher conserved charges of single solitons by Freeman [-@freeman1995]. Complex-energy and negative-energy solutions are in Khastgir and Sasaki [-@khastgir1996], and the complex time delays and singular final states of $a_2^{(1)}$ in Takács and Watts [-@tw1999]. The missing charges were found by McGhee [-@mcghee1994]. The formulation as a selection rule from $\Theta$ is new here.
  ```
  with:
  ```markdown
  Reality of soliton energies was shown by Olive, Turok and Underwood [-@olive1993a; -@olive1993b] and Zhu and Caldi [@zhu1995], and of the higher conserved charges of single solitons by Freeman [-@freeman1995]. Using antilinear $\mathcal{PT}$ symmetry together with asymptotic additivity to establish reality of conserved charges for multi-soliton solutions at real rapidities was introduced by Cen and Fring [@cen2016] and Cen, Correa and Fring [@cen2017], and applied to Bullough–Dodd solitons and breathers by Correa, Fring and Taira [@correa2022] (see also [@fring2020]). Complex-energy and negative-energy solutions are in Khastgir and Sasaki [-@khastgir1996], and the complex time delays and singular final states of $a_2^{(1)}$ in Takács and Watts [-@tw1999]. The missing charges were found by McGhee [-@mcghee1994]. The formulation of @prp-theta-selection-rule as a general selection rule on multicomponent Hirota data $\{(a_i,\theta_i)\}$ under $\Theta:(a,\theta)\mapsto(\bar a,\bar\theta)$, including complex-rapidity conjugate pairs and the charge-conjugation sign $\epsilon_Q=(-1)^{s-1}$ for higher charges, is new here.
  ```
- In `appendices/E-original-results.qmd`, line 8, replace:
  ```markdown
  - **[Theorem, given additivity of the higher charges]** A selection rule from the antilinear symmetry $\Theta$: a Hirota solution whose asymptotic data are $\Theta$-invariant has $Q=\epsilon_Q\bar Q$ for every conserved charge, and real energy and momentum whenever the multiset of $(M_{a_i},\theta_i)$ is self-conjugate (@sec-reality-energy-momentum).
  ```
  with:
  ```markdown
  - **[Theorem, given additivity of the higher charges]** A selection rule on multicomponent Hirota data from the antilinear symmetry $\Theta$, extending the $\mathcal{PT}$ and asymptotic-additivity mechanism of [@cen2016; @cen2017] to complex-rapidity conjugate pairs and $C$-odd higher charges: a regular Hirota multi-soliton solution whose constituent data $\{(a_i,\theta_i)\}$ are $\Theta$-invariant under $(a,\theta)\mapsto(\bar a,\bar\theta)$ has $Q=\epsilon_Q\bar Q$ for every additive conserved charge, and real positive energy and real momentum whenever the multiset of $(M_{a_i},\theta_i)$ is self-conjugate with $\operatorname{Re}\cosh\theta_i>0$ (@prp-theta-selection-rule).
  ```

---

## Item 2. `@prp-hirota-energy` (`11-reality.qmd:32–41`)

### Overlap

- **Overlap with Phase 1:** **in both**.
  - **Covered in Phase 1:** Rank 1 (`11-reality.qmd:32–41`) and Rank 4 (`10-solitons.qmd:174`). Phase 1 read `olive1993b` (`hep-th/9305160`), `zhu1995` (`hep-th/9404097`), and `tw1999` (`hep-th/9810006`), checked the pointwise identity $T_{\mu\nu}^{\rm can} = (\eta_{\mu\nu}\partial^2 - \partial_\mu\partial_\nu)C$, derived the boundary conditions $\operatorname{Re}\cosh\theta_i > 0$ and $\prod_{k<l}A_{a_k a_l}(\theta_k-\theta_l)\notin\{0,\infty\}$ required for the surface term $[\partial_x\ln\tau_j]_{-\infty}^{+\infty}$ to evaluate to $\sum_i M_{a_i}\cosh\theta_i$, and demonstrated failure when $\operatorname{Re}\cosh\theta_i < 0$ or $A_{12}(2iu) = 0$.
  - **New in Phase 2:** Direct numerical integration of the canonical energy and momentum densities of the actual field $\vec\phi$ on an $a_2^{(1)}$ breather (`check_phase2_items.py` Item 2).

### Check

1. **Exact source check of `[@olive1993a; @olive1993b]`, `[@zhu1995]`, and `[@tw1999, eq. 4.6]`:**
   - **Olive, Turok and Underwood (`hep-th/9305160`, `olive1993b`):**
     - In Section 8.2 (eqs. (8.6)–(8.7), lines 1295–1325), they quote from `olive1993a` (*Nucl. Phys.* B401 (1993) 663) that on any soliton solution the energy–momentum tensor is improved to a total derivative with $P^\pm = \pm\frac{\sqrt{2}}{\beta^2}\big[(\partial_+ - \partial_-)\sum_{i=0}^r \frac{2}{\alpha_i^2}n_i\ln\tau_i\big]_{-\infty}^\infty$ (in their convention where $\tau_i$ is the $n_i$-th power of the level-1 $\tau_i$).
     - In Section 8.3 ("Continuously varying soliton masses?", lines 1330–1354), they examine complex rapidities ($\sigma_1 = \mu e^{\eta_1}$ complex) and write (lines 1335–1348):
       > *"Notice that in the single soliton mass calculation using (8.7) only the dominant behaviour of $\ln\tau_i$ as $x$ approaches $\pm\infty$ enters. In fact this argument holds good even if $e^{\eta_i}$ is complex, **providing that its real part is positive** ... Using our general formula (3.11) we write the two-soliton solution depending on $e^{\eta_1}$ and $e^{\eta_2}$. We take these to be complex conjugate variables $e^{\eta_1} = (e^{\eta_2})^* = e^{\eta + i\theta}$, where $\eta$ and $\theta$ are real and $-\pi/2 < \theta < \pi/2$, in order to satisfy the reality condition discussed above."*
   - **Takács and Watts (`hep-th/9810006`, `tw1999`):**
     - In Section 4.1 (lines 631–636, eq. (4.6)), they write for $a_2^{(1)}$:
       > *"The energy and momentum of a multi-soliton solution with rapidities $\theta_i$ turns out to be given completely generally by the formula [30]: $E = M\sum_i\cosh\theta_i,\quad P = M\sum_i\sinh\theta_i,\quad M = \frac{6\sqrt{3}m}{\beta^2}$ (4.6) where the $\theta_i$ are the (possibly complex) rapidities of the constituent solitons."*
     - Then at lines 686–698 (eqs. (4.10)–(4.11)), they sketch the derivation from $\tau_j = 1 + \dots + (\prod_{m<n}A^{(a_m a_n)})\prod_m W_m$, which explicitly assumes that as $x\to+\infty$ every $W_m \propto \exp(\sqrt{3}\cosh\theta_m\cdot x)$ grows ($\operatorname{Re}\cosh\theta_m > 0$) and that the top coefficient $\prod_{m<n}A^{(a_m a_n)}(\theta_m-\theta_n)$ is non-zero and finite.
   - **Zhu and Caldi (`hep-th/9404097`, `zhu1995`):**
     - In Section 3 (lines 245–368), Zhu and Caldi compute the **local energy density** $\mathcal H(x,t)$ of single solitons and breathers in $a_n^{(1)}$ and study when the density itself is pointwise real and positive (which occurs for static 1-solitons only at discrete phase values). They do not prove the general complex-rapidity multi-soliton sum formula (`eq. 4.6` of `tw1999`).

2. **Derivation of $E = \sum_i M_{a_i}\cosh\theta_i$ and $P = \sum_i M_{a_i}\sinh\theta_i$ from `eq-hirota-bilinear` and `C`, and where complex rapidities fail:**
   - With $C = -\frac{2}{\beta^2}\sum_{j=0}^r \frac{2}{\vec\alpha_j^2}\ln\tau_j$, we verified in Phase 1 Rank 4 (`check_ch11_reality_breathers.py` Check 1) that on every Hirota multi-soliton solution the canonical stress-energy tensor $T_{\mu\nu} = \partial_\mu\vec\phi\cdot\partial_\nu\vec\phi - \eta_{\mu\nu}\mathcal L$ satisfies
     $$T_{00} = -\partial_x^2 C = \frac{2}{\beta^2}\sum_{j=0}^r \frac{2}{\vec\alpha_j^2}\partial_x^2\ln\tau_j,\qquad T_{01} = -\partial_t\partial_x C = \frac{2}{\beta^2}\sum_{j=0}^r \frac{2}{\vec\alpha_j^2}\partial_t\partial_x\ln\tau_j$$
     pointwise.
   - Provided no $\tau_j(x,t)$ vanishes on the integration line $x\in(-\infty,\infty)$, integrating over $x$ gives the exact surface terms:
     $$E = \int_{-\infty}^\infty T_{00}\,dx = \frac{2}{\beta^2}\sum_{j=0}^r \frac{2}{\vec\alpha_j^2}\big[\partial_x\ln\tau_j\big]_{-\infty}^{+\infty},\qquad P = -\int_{-\infty}^\infty T_{01}\,dx = -\frac{2}{\beta^2}\sum_{j=0}^r \frac{2}{\vec\alpha_j^2}\big[\partial_t\ln\tau_j\big]_{-\infty}^{+\infty}.$$
   - For an $N$-soliton solution (`eq-n-soliton`), $\tau_j = \sum_{S\subseteq\{1,\dots,N\}} \big(\prod_{k\in S}\omega^{j a_k}E_k\big)\prod_{k<l\in S}A_{a_k a_l}(\theta_k-\theta_l)$ with $E_k = \exp\big(m_{a_k}(x\cosh\theta_k - t\sinh\theta_k) + \xi_k\big)$.
     - **Condition 1 ($\operatorname{Re}\cosh\theta_k > 0$):** As $x\to -\infty$, $E_k \to 0$ for all $k$ if and only if $\operatorname{Re}\cosh\theta_k > 0$ for all $k$ (i.e. $\operatorname{Im}\theta_k \in (-\frac\pi2,\frac\pi2)\pmod{2\pi}$), in which case $\tau_j \to 1$ and $\partial_\mu\ln\tau_j \to 0$ at $x\to -\infty$. If a constituent has $\operatorname{Re}\cosh\theta_k < 0$, then $E_k\to\infty$ as $x\to -\infty$ and $E_k\to 0$ as $x\to +\infty$, which **flips the sign** of its contribution to $+\operatorname{sgn}(\operatorname{Re}\cosh\theta_k)\,M_{a_k}\cosh\theta_k$ (so $\operatorname{Re}E > 0$ even when $\operatorname{Re}\cosh\theta_k < 0$, contrary to simply plugging $\operatorname{Im}\theta_k = \pi$ into $M_a\cosh\theta$).
     - **Condition 2 ($\prod_{k<l}A_{a_k a_l}(\theta_k-\theta_l)\notin\{0,\infty\}$):** When $\operatorname{Re}\cosh\theta_k > 0$ for all $k$ AND the full $N$-soliton coefficient $A_{\{1,\dots,N\}} = \prod_{1\le k<l\le N}A_{a_k a_l}(\theta_k-\theta_l)$ is non-zero and finite, the single term $S = \{1,\dots,N\}$ dominates as $x\to +\infty$:
       $$\tau_j(x,t) \underset{x\to+\infty}{\sim} A_{\{1,\dots,N\}}\prod_{k=1}^N \omega^{j a_k}E_k \implies \partial_x\ln\tau_j \underset{x\to+\infty}{\longrightarrow} \sum_{k=1}^N m_{a_k}\cosh\theta_k,\quad -\partial_t\ln\tau_j \underset{x\to+\infty}{\longrightarrow} \sum_{k=1}^N m_{a_k}\sinh\theta_k,$$
       which yields $E = \sum_{k=1}^N M_{a_k}\cosh\theta_k$ and $P = \sum_{k=1}^N M_{a_k}\sinh\theta_k$.
       If instead $A_{\{1,\dots,N\}} = 0$ — which happens for an $a_n^{(1)}$ breather $(a, h-a)$ at $u = \frac\pi2 - \frac{\pi a}{h} \in (0,\frac\pi2)$ whenever $a \ne h/2$ (e.g. $u = \pi/6$ in $a_2^{(1)}$) — the $E_1 E_2$ term in $\tau_j$ vanishes identically, and $\tau_j = 1 + \omega^{ja}E_1 + \omega^{-ja}E_2$ grows as $E_1 + E_2 \sim e^{m_a x\cosh\theta_0\cos u}$, giving **half** the breather energy ($M_a\cosh\theta_0\cos u$ instead of $2M_a\cosh\theta_0\cos u$).

3. **Numerical integration of the canonical energy and momentum densities of $\vec\phi$ on an $a_2^{(1)}$ breather (`check_phase2_items.py` Item 2):**
   We evaluated the canonical energy density $T_{00} = \frac12(\partial_t\vec\phi)^2 + \frac12(\partial_x\vec\phi)^2 + V_{\rm imag}(\vec\phi)$ and momentum density $\mathcal P = -\partial_t\vec\phi\cdot\partial_x\vec\phi$ directly from the complex field $\vec\phi(x,t)$ of an $a_2^{(1)}$ breather at $\theta_0 = 0.3$, $u = 1.0 \in (\pi/6,\pi/2)$, $\xi_1 = 0.2 + 0.4i$, $\xi_2 = \bar\xi_1 + i\pi$, $t = 0.2$ ($m=\beta=1, M_1 = 6\sqrt{3}$):
   ```text
   a_2^(1) breather (th0=0.3, u=1.0):
     Canonical int T00 dx = 11.7391228112 + -1.70e-11j,  2 M1 cosh(th0) cos(u) = 11.7391228113, diff = 1.21e-10
     Canonical int P   dx = 3.4197545339 + -1.98e-10j,  2 M1 sinh(th0) cos(u) = 3.4197545341, diff = 1.70e-10
   ```

### Verdict

- **Verdict:** **wrong** (`@prp-hirota-energy` omits the hypotheses $\operatorname{Re}\cosh\theta_i > 0$, $\prod_{k<l}A_{a_k a_l}(\theta_k-\theta_l)\notin\{0,\infty\}$, and $\tau_j\ne 0$ on $\mathbb R$, which are explicitly required by the derivation and stated in `olive1993b`, Section 8.3).

### Paste-ready replacement text

- In `part3-classical/11-reality.qmd`, lines 32–41, replace:
  ```markdown
  ::: {#prp-hirota-energy}
  **[Established; reformulated]** For a Hirota solution, the energy and momentum densities are total derivatives [@olive1993a; @olive1993b]. They integrate to

  $$
  E=\sum_iM_{a_i}\cosh\theta_i,\qquad P=\sum_iM_{a_i}\sinh\theta_i,
  $$

  where $(a_i,\theta_i)$ are the species and (possibly complex) rapidities of the constituent solitons [@tw1999, eq. 4.6]. The masses $M_a$ are real and positive. For the higher conserved charges of a single soliton the same holds [@freeman1995]; for multi-soliton solutions it is assumed here, on the basis of additivity (see the remark below).

  :::
  ```
  with:
  ```markdown
  ::: {#prp-hirota-energy}
  **[Established; reformulated]** For a Hirota multi-soliton solution whose constituent solitons $(a_i,\theta_i)$ satisfy $\operatorname{Re}\cosh\theta_i>0$ and $\prod_{k<l}A_{a_ka_l}(\theta_k-\theta_l)\notin\{0,\infty\}$, and for which no $\tau_j$ vanishes on the real line, the energy and momentum densities are total derivatives [@olive1993a; @olive1993b] that integrate to

  $$
  E=\sum_iM_{a_i}\cosh\theta_i,\qquad P=\sum_iM_{a_i}\sinh\theta_i
  $$

  [@olive1993b, sec. 8.3; @tw1999, eq. 4.6]. The masses $M_a$ are real and positive. (If $\operatorname{Re}\cosh\theta_i<0$, the limits $x\to\pm\infty$ exchange roles for $E_i$ and its contribution acquires a minus sign, so $\operatorname{Re}E$ remains positive.) For the higher conserved charges of a single soliton the same holds [@freeman1995]; for multi-soliton solutions it is assumed here, on the basis of additivity (see the remark below).
  :::
  ```

---

## Item 3. The additivity assumption tested (`11-reality.qmd:46`)

### Overlap

- **Overlap with Phase 1:** **in both**.
  - **Covered in Phase 1:** Rank 2 (`11-reality.qmd:46`). Phase 1 read Freeman (`hep-th/9408092`), constructed the $a_2^{(1)}$ spin-2 current symbolically in `referee_scratch/classical_09_11/check_ch11_reality_breathers.py` (Check 4), evaluated it on single solitons of species $a=1,2$ to confirm $q_2 = -q_1$ (eigenvector of the $A_2$ Cartan matrix with eigenvalue $3$), and integrated it numerically on an $a_2^{(1)}$ breather at complex conjugate rapidities $\theta_0\pm iu$ at two different times $t=0.0$ and $t=0.7$, confirming $Q_2 = 2iq_1 e^{2\theta_0}\sin 2u$ to $5\times 10^{-14}$!

### Check

1. **Exact verification of Freeman (`hep-th/9408092`, `freeman1995`):**
   - In `hep-th/9408092`, Section 3.1, eq. (3.20) (lines 475–485), Freeman proves that on single- and multi-soliton solutions, the conserved charge density for spin $N$ is a total derivative:
     $$Q_{-N} = m^N\left[\frac{\langle\Lambda_0|g\hat E_{-N}|\Lambda_0\rangle}{\langle\Lambda_0|g|\Lambda_0\rangle}\right]_{x^-=-\infty}^\infty.$$
   - In Section 3.2, eqs. (3.44), (3.47)–(3.49) (lines 678–709), Freeman evaluates this surface term for single solitons of an arbitrary untwisted affine Toda theory $\hat g$ and concludes:
     > *"Thus we arrive finally at the result that the values of the conserved charge $Q_{-N}$ for soliton solutions form a left eigenvector of the Cartan matrix of $g$, with eigenvalue $2(1-\cos N\pi/h)$. By a suitable choice of phase all of these values can be taken to be real. The antisoliton solutions can be treated similarly; these are related to the solitons by $Q_{-N}(\text{antisoliton}) = (-1)^{N-1}Q_{-N}(\text{soliton})$ (3.49)."*
   - For $a_2^{(1)}$ ($h=3$) and spin $N=2$, the eigenvalue of the $A_2$ Cartan matrix $\begin{pmatrix}2&-1\\-1&2\end{pmatrix}$ is $2(1-\cos(2\pi/3)) = 3$, whose eigenvector is $(q_1, q_2) = q_1(1, -1)$, so $q_2 = -q_1$, exactly as stated in `11-reality.qmd:46`.

2. **Direct numerical test of additivity on an $a_2^{(1)}$ breather (`referee_scratch/classical_09_11/check_ch11_reality_breathers.py`, Check 4):**
   - Using the spin-2 conserved density $J_{Q_2} = J_{+3} + J_{-1}$ (with $\partial_+\leftrightarrow\partial_-$ so that $Q_2$ scales as $e^{+2\theta}$) constructed in Item 1 above, we integrated $J_{Q_2}(\vec\phi(x,t))$ numerically over $x\in[-30,30]$ on single solitons ($a=1,2$) and on the regular $a_2^{(1)}$ breather $\{(1, \theta_0+iu), (2, \theta_0-iu)\}$ with $\theta_0 = 0.3, u = 1.0 \in (\pi/6,\pi/2), \xi_1 = 0.1+0.2i, \xi_2 = \bar\xi_1 + i\pi$ at $t=0.0$ and $t=0.7$:
     ```text
     th=0.0 [d_- (scales as e^{+2 th})]: Q(a=1) = 38.183766-0.000000j, Q(a=2) = -38.183766+0.000000j
     th=0.4 [d_- (scales as e^{+2 th})]: Q(a=1) = 84.979534-0.000000j, Q(a=2) = -84.979534+0.000000j, Q(a=1)/e^{2th} = 38.183766
     Breather (th0=0.3, u=1.0, t=0.0): Q_2(num) = 0.000000+126.529388j, 2i q_1 e^{2 th0} sin(2u) = 0.000000+126.529388j, diff = 5.00e-14
     Breather (th0=0.3, u=1.0, t=0.7): Q_2(num) = 0.000000+126.529388j, 2i q_1 e^{2 th0} sin(2u) = 0.000000+126.529388j, diff = 6.81e-14
     ```
   - **Analytical reason why additivity holds at complex rapidities:** Freeman's surface-term formula (3.20), $Q_{-N} = m^N\big[\langle\Lambda_0|g\hat E_{-N}|\Lambda_0\rangle / \langle\Lambda_0|g|\Lambda_0\rangle\big]_{-\infty}^\infty$, holds for any multi-soliton group element $g = \prod_{k=1}^n \exp(Q_k \hat F(a_k, \rho_k))$. Because $[\hat E_{-N}, \hat F(a_k, \rho_k)] = (\alpha_{a_k}, E_{h-\nu})\rho_k^{-N}\hat F(a_k, \rho_k)$ acts as a derivation on the product $\prod_k \exp(Q_k\hat F_k)$, whenever $\operatorname{Re}\cosh\theta_k > 0$ for all $k$ and $\prod_{k<l}A_{a_k a_l}(\theta_k-\theta_l)\notin\{0,\infty\}$, the top term where all $\hat F_k$ act to their nilpotency bound $n_0=1$ dominates as $x\to +\infty$, so the surface term is strictly additive: $Q_{-N} = \sum_{k=1}^n Q_{-N}^{(1)}(a_k,\theta_k)$!

### Verdict

- **Verdict:** **correct** (the $a_2^{(1)}$ spin-2 calculation, Freeman citation, and breather formula $2iq_1 e^{2\theta_0}\sin 2u$ at `11-reality.qmd:46` are completely accurate; moreover, our numerical check and Freeman's eq. (3.20) confirm that additivity actually holds on regular breathers).

---

## Item 4. Time delays and the Wigner sign (`10-solitons.qmd:122`)

### Overlap

- **Overlap with Phase 1:** **in both**.
  - **Covered in Phase 1:** Rank 13 (`10-solitons.qmd:122`) and Section 4 Observation 1 (overload of $E_1$). Phase 1 checked `fring1994` (`hep-th/9405034`), derived the asymptotic trajectories of soliton 1 as $t\to\pm\infty$, verified $\Delta x > 0$, $\Delta t = -\Delta x/v_1 < 0$, $E_1\Delta x = -p_1\Delta t = -(2h/\beta^2)\ln A_{a_1 a_2}$, Wigner's relation $d\delta_{ab}/d\theta = p_1\Delta t = (2h/\beta^2)\ln A_{ab}(\theta) < 0$, and the sine-Gordon limit $(16/\beta_{\rm SG}^2)\ln\tanh\frac\theta2$ from `eq-sg-s0`.
  - **New in Phase 2:**
    1. Explicit verification of the trajectory shift for **soliton 2** as well as soliton 1 (confirming the caption of `@fig-two-soliton-shift`, `10-solitons.qmd:124`, that *each* soliton is displaced forward along its own direction of motion).
    2. Comparison of $(2h/\beta^2)\ln A_{ab}(\theta)$ against the $\beta\to 0$ limit of Hollowood's $a_n^{(1)}$ quantum soliton S-matrix scalar factor from `@sec-an-example` (`21-construction.qmd:43–55`) for $a_2^{(1)}$ and $a_3^{(1)}$ (`check_phase2_items.py` Item 4).

### Check

1. **Asymptotic trajectories of both soliton 1 and soliton 2 (`10-solitons.qmd:122, 124`):**
   - Consider the two-soliton $\tau$-function (`eq-n-soliton` for $N=2$):
     $$\tau_j(x,t) = 1 + \omega^{j a_1}E_1 + \omega^{j a_2}E_2 + \omega^{j(a_1+a_2)}A_{a_1 a_2}(\theta_1-\theta_2)E_1 E_2,$$
     with $E_k = \exp\big(m_{a_k}(x\cosh\theta_k - t\sinh\theta_k) + \xi_k\big) = \exp\big(m_{a_k}\cosh\theta_k(x - v_k t) + \xi_k\big)$, where $v_k = \tanh\theta_k$ and $\theta_1 > \theta_2$ ($v_1 > v_2$).
   -For real rapidities $\theta_1 \ne \theta_2$, $\cosh(\theta_1-\theta_2) > 1$ and $\cos\frac{\pi(a_1-a_2)}{h} > \cos\frac{\pi(a_1+a_2)}{h}$ for $1\le a_1,a_2\le h-1$, so
     $$0 < A_{a_1 a_2}(\theta_1-\theta_2) = \frac{\cosh(\theta_1-\theta_2)-\cos\frac{\pi(a_1-a_2)}{h}}{\cosh(\theta_1-\theta_2)-\cos\frac{\pi(a_1+a_2)}{h}} < 1 \implies \ln A_{a_1 a_2}(\theta_1-\theta_2) < 0.$$
   - **Follow soliton 1 ($x = v_1 t + x_0$):**
     Along this worldline, $E_2 \propto \exp\big(m_{a_2}\cosh\theta_2(v_1 - v_2)t\big)$. Since $v_1 - v_2 > 0$:
     - As $t\to -\infty$: $E_2 \to 0$, so $\tau_j \simeq 1 + \omega^{j a_1}E_1$, centred at $x_1^{\rm in}(t) = v_1 t - \frac{\operatorname{Re}\xi_1}{m_{a_1}\cosh\theta_1}$.
     - As $t\to +\infty$: $E_2 \to \infty$, so $\tau_j \simeq \omega^{j a_2}E_2\big(1 + \omega^{j a_1}A_{a_1 a_2}E_1\big)$, centred at $x_1^{\rm out}(t) = v_1 t - \frac{\operatorname{Re}\xi_1 + \ln A_{a_1 a_2}}{m_{a_1}\cosh\theta_1}$.
     - Thus soliton 1 is shifted in $x$ by
       $$\Delta x_1 = x_1^{\rm out}(t) - x_1^{\rm in}(t) = -\frac{\ln A_{a_1 a_2}(\theta_1-\theta_2)}{m_{a_1}\cosh\theta_1} > 0.$$
       When $\theta_1 > 0$ ($v_1 > 0$, moving right), $\Delta x_1 > 0$ is forward along its direction of motion, and the time shift at a fixed position $x$ is $\Delta t_1 = t_1^{\rm out}(x) - t_1^{\rm in}(x) = -\Delta x_1 / v_1 < 0$ (a time advance).
   - **Follow soliton 2 ($x = v_2 t + x_0$):**
     Along this worldline, $E_1 \propto \exp\big(m_{a_1}\cosh\theta_1(v_2 - v_1)t\big)$. Since $v_2 - v_1 < 0$:
     - As $t\to -\infty$: $E_1 \to +\infty$, so $\tau_j \simeq \omega^{j a_1}E_1\big(1 + \omega^{j a_2}A_{a_1 a_2}E_2\big)$, centred at $x_2^{\rm in}(t) = v_2 t - \frac{\operatorname{Re}\xi_2 + \ln A_{a_1 a_2}}{m_{a_2}\cosh\theta_2}$.
     - As $t\to +\infty$: $E_1 \to 0$, so $\tau_j \simeq 1 + \omega^{j a_2}E_2$, centred at $x_2^{\rm out}(t) = v_2 t - \frac{\operatorname{Re}\xi_2}{m_{a_2}\cosh\theta_2}$.
     - Thus soliton 2's coordinate shift is
       $$\Delta x_2 = x_2^{\rm out}(t) - x_2^{\rm in}(t) = +\frac{\ln A_{a_1 a_2}(\theta_1-\theta_2)}{m_{a_2}\cosh\theta_2} < 0.$$
       In `@fig-two-soliton-shift` (`10-solitons.qmd:124`), $\theta_1 = 0.6 > 0$ (moving right, $+x$) and $\theta_2 = -0.6 < 0$ (moving left, $-x$). Since soliton 2 moves in the $-x$ direction, a shift $\Delta x_2 < 0$ is **forward along its own direction of motion** by the positive distance $\lvert\Delta x_2\rvert = -\frac{\ln A_{12}}{m_{a_2}\cosh\theta_2} > 0$, and its time shift is $\Delta t_2 = -\Delta x_2 / v_2 < 0$ (also a time advance). This confirms both `10-solitons.qmd:122` and the caption of `@fig-two-soliton-shift` (`10-solitons.qmd:124`), as well as Fring, Johnson, Kneipp and Olive (`hep-th/9405034`, eqs. (2.8)–(2.14)).

2. **Wigner's relation and comparison with sine-Gordon (`03-sine-gordon.qmd:136`) and Hollowood's $a_n^{(1)}$ S-matrix (`21-construction.qmd:43–55`):**
   - Using $M_{a_1} = \frac{2h}{\beta^2}m_{a_1}$ (`eq-soliton-mass`), soliton 1's energy $\mathcal E_1 = M_{a_1}\cosh\theta_1$ and momentum $p_1 = M_{a_1}\sinh\theta_1$ satisfy
     $$\mathcal E_1\Delta x_1 = -p_1\Delta t_1 = -\frac{2h}{\beta^2}\ln A_{a_1 a_2}(\theta_1-\theta_2) > 0.$$
   - In the rest frame of soliton 2 ($\theta_2 = 0$, $\theta = \theta_1 > 0$), $d\mathcal E_1 = M_{a_1}\sinh\theta\,d\theta = p_1\,d\theta$, so Wigner's time-delay formula $\Delta t_1 = \frac{d\delta_{a_1 a_2}}{d\mathcal E_1} = \frac{1}{p_1}\frac{d\delta_{a_1 a_2}}{d\theta}$ gives
     $$\frac{d\delta_{ab}(\theta)}{d\theta} = p_1\Delta t_1 = \frac{2h}{\beta^2}\ln A_{ab}(\theta) < 0.$$
   - **Sine-Gordon check (`eq-sg-s0`, `03-sine-gordon.qmd:136`):**
     For $a_1^{(1)}$ ($h=2$, $a=b=1$, $\beta_{\rm SG}^2 = 2\beta^2$), $A_{11}(\theta) = \frac{\cosh\theta - 1}{\cosh\theta + 1} = \tanh^2\frac{\theta}{2}$, so
     $$\frac{2h}{\beta^2}\ln A_{11}(\theta) = \frac{4}{\beta^2}\ln\tanh^2\frac{\theta}{2} = \frac{16}{\beta_{\rm SG}^2}\ln\tanh\frac{\theta}{2}.$$
     From `eq-sg-s0`, $-i\ln(-S_0(\theta)) = -\int_0^\infty \frac{dt}{t}\frac{\sinh\frac{(\pi-\xi)t}{2}\sin\theta t}{\sinh\frac{\xi t}{2}\cosh\frac{\pi t}{2}}$. Differentiating twice with respect to $\theta$ and taking $\xi \simeq \beta_{\rm SG}^2/8 \to 0$ gives
     $$\frac{d}{d\theta}\left(-i\frac{d\ln S_0}{d\theta}\right) \underset{\xi\to 0}{\simeq} \frac{2}{\xi}\int_0^\infty dt\,\tanh\frac{\pi t}{2}\,\sin\theta t = \frac{2}{\xi}\frac{1}{\sinh\theta} = \frac{16}{\beta_{\rm SG}^2}\frac{d}{d\theta}\ln\tanh\frac{\theta}{2},$$
     which integrates (with vanishing boundary condition at $\theta\to+\infty$) to $-i\frac{d\ln S_0}{d\theta} \simeq \frac{16}{\beta_{\rm SG}^2}\ln\tanh\frac{\theta}{2}$.
   - **Comparison with Hollowood's $a_n^{(1)}$ quantum soliton S-matrix (`21-construction.qmd:43–55`, `hollowood1993b`):**
     In `@sec-an-example` (`21-construction.qmd:43–55`), the scalar factor $S_{ab}(\theta)$ of Hollowood's $U_q(a_n^{(1)})$ soliton S-matrix has the integral representation (Fring–Johnson–Kneipp–Olive `hep-th/9405034`, Section 4) whose $\beta\to 0$ limit gives
     $$\lim_{\beta\to 0}\frac{\beta^2}{2h}\frac{d}{d\theta}\left(-i\frac{d\ln S_{ab}(\theta)}{d\theta}\right) = 2\int_0^\infty dt\,\frac{\cosh((\pi-\phi_-)t) - \cosh((\pi-\phi_+)t)}{\sinh(\pi t)}\sin\theta t = \frac{d}{d\theta}\ln A_{ab}(\theta),$$
     where $\phi_- = \pi\lvert a-b\rvert/h$ and $\phi_+ = \pi(a+b)/h$. Verified numerically in `check_phase2_items.py` (Item 4) for all $(a,b)$ in $a_2^{(1)}$ and $a_3^{(1)}$ at $\theta = 1.1$:
     ```text
     a_2^(1) (a=1,b=1): d/dth ln A_ab = 1.3819950821, semiclassical S-matrix integral = 1.3819950821, diff = 2.22e-16
     a_2^(1) (a=2,b=1): d/dth ln A_ab = 0.6425061154, semiclassical S-matrix integral = 0.6425061154, diff = 0.00e+00
     a_3^(1) (a=1,b=1): d/dth ln A_ab = 1.1974222962, semiclassical S-matrix integral = 1.1974222962, diff = 0.00e+00
     a_3^(1) (a=2,b=1): d/dth ln A_ab = 0.8270265648, semiclassical S-matrix integral = 0.8270265648, diff = 2.22e-16
     a_3^(1) (a=2,b=2): d/dth ln A_ab = 1.4974011068, semiclassical S-matrix integral = 1.4974011068, diff = 4.44e-16
     a_3^(1) (a=3,b=1): d/dth ln A_ab = 0.2999788106, semiclassical S-matrix integral = 0.2999788106, diff = 5.55e-17
     ```

### Verdict

- **Verdict:** **correct** (all signs, Wigner's relation, the trajectory shifts of both solitons, and the sine-Gordon / Hollowood semiclassical limits are completely accurate; see Section 6 Observation 1 for the minor notation collision where $E_1$ denotes both the exponential of soliton 1 and the energy of soliton 1 in `10-solitons.qmd:122`).

---

## Item 5. Lax pair, zero curvature, boost invariance and conserved-charge spins (`09-lagrangian.qmd:32–99`)

### Overlap

- **Overlap with Phase 1:** **in both**.
  - **Covered in Phase 1:** Rank 8 (`09-lagrangian.qmd:90–96`, `eq-spin3-current`, `exr-spin3`), Rank 10 (`09-lagrangian.qmd:84–88`, exponents modulo $kh$), and Rank 16 (`09-lagrangian.qmd:57–79`, `eq-lax-pair`, `eq-zero-curvature`, boost identity, `exr-lax`).
  - **New in Phase 2:**
    1. Explicit $3\times 3$ matrix check of `eq-lax-pair` and `eq-zero-curvature` for the non-simply-laced twisted algebra $a_2^{(2)}$ (Bullough–Dodd) in `check_phase2_items.py` (Item 5).
    2. Explicit derivation showing why a boost requires a constant gauge transformation in a non-principal gradation $(s_0,\dots,s_r)\ne(1,\dots,1)$.
    3. Explicit check of the centralizer of $E_1$ in $a_4^{(2)}$ (`check_phase2_items.py` Item 5) in addition to $a_2^{(2)}, d_4^{(3)}, e_6^{(2)}, c_2^{(1)}, d_3^{(2)}, a_n^{(1)}$ checked in `check_rank10_exponents.py`.
    4. Literature check of `[@feigin1993; @dgz1992; @niedermaier1994]` at `09-lagrangian.qmd:32`.

### Check

1. **Zero-curvature calculation (`eq-lax-pair`, `eq-zero-curvature`, `09-lagrangian.qmd:59–73`) for general $\hat g$ and explicit $3\times 3$ check for $a_2^{(2)}$:**
   - Using $[E_{\vec\alpha_i}, E_{-\vec\alpha_j}] = \delta_{ij}\frac{2\vec\alpha_j}{\vec\alpha_j^2}\cdot\vec H$ and $[\vec H, E_{\pm\vec\alpha_j}] = \pm\vec\alpha_j E_{\pm\vec\alpha_j}$, we have:
     $$\partial_+\mathcal A_- - \partial_-\mathcal A_+ = -\beta\,\partial_+\partial_-\vec\phi\cdot\vec H - \frac{m\beta}{2\lambda}\sum_j c_j(\vec\alpha_j\cdot\partial_+\vec\phi)e^{\beta\vec\alpha_j\cdot\vec\phi/2}E_{-\vec\alpha_j} - \frac{m\beta\lambda}{2}\sum_j c_j(\vec\alpha_j\cdot\partial_-\vec\phi)e^{\beta\vec\alpha_j\cdot\vec\phi/2}E_{\vec\alpha_j},$$
     $$[\mathcal A_+, \mathcal A_-] = \frac{m\beta}{2\lambda}\sum_j c_j(\vec\alpha_j\cdot\partial_+\vec\phi)e^{\beta\vec\alpha_j\cdot\vec\phi/2}E_{-\vec\alpha_j} + \frac{m\beta\lambda}{2}\sum_j c_j(\vec\alpha_j\cdot\partial_-\vec\phi)e^{\beta\vec\alpha_j\cdot\vec\phi/2}E_{\vec\alpha_j} - m^2\sum_j c_j^2 e^{\beta\vec\alpha_j\cdot\vec\phi}\frac{2\vec\alpha_j}{\vec\alpha_j^2}\cdot\vec H.$$
     Adding the two equations, the $E_{\pm\vec\alpha_j}$ terms cancel identically, and with $c_j^2 = n_j\vec\alpha_j^2/2$ the Cartan coefficient is
     $$-\left(\beta\,\partial_+\partial_-\vec\phi + m^2\sum_{j=0}^r n_j\vec\alpha_j e^{\beta\vec\alpha_j\cdot\vec\phi}\right)\cdot\vec H = 0,$$
     which is equivalent to `eq-toda-eom` for every untwisted and twisted affine algebra $\hat g$.
   - **Explicit $3\times 3$ matrix check for $a_2^{(2)}$ (Bullough–Dodd) (`check_phase2_items.py` Item 5):**
     With $\alpha_1 = 1/\sqrt{2}$ ($n_1=2, c_1 = 1/\sqrt{2}$), $\alpha_0 = -\sqrt{2}$ ($n_0=1, c_0 = 1$), $H = \frac{1}{\sqrt{2}}\operatorname{diag}(1,0,-1)$, $e_1 = E_{12}+E_{23}$, $f_1 = 2(E_{21}+E_{32})$, $e_0 = E_{31}$, $f_0 = E_{13}$, SymPy confirms that $\partial_+\mathcal A_- - \partial_-\mathcal A_+ + [\mathcal A_+,\mathcal A_-] = 0_{3\times 3}$ if and only if $\partial_+\partial_-\phi = -\frac{m^2}{\beta}\big(\sqrt{2}e^{\beta\phi/\sqrt{2}} - \sqrt{2}e^{-\sqrt{2}\beta\phi}\big)$:
     ```text
     a_2^(2) (Bullough-Dodd) 3x3 zero-curvature residual with EOM substituted: True
     ```

2. **Boost transformation of $\mathcal A_\pm$ (`09-lagrangian.qmd:77`) and why another gradation requires a constant gauge transformation:**
   - Under $\vec\phi_\vartheta(x^+,x^-) = \vec\phi(e^{-\vartheta}x^+, e^{\vartheta}x^-)$, the chain rule gives $(\partial_\pm\vec\phi_\vartheta)(x) = e^{\mp\vartheta}(\partial_\pm\vec\phi)(x')$ where $x' = (e^{-\vartheta}x^+, e^\vartheta x^-)$.
   - In the principal gradation, $\lambda$ multiplies every $E_{\vec\alpha_j}$ ($j=0,\dots,r$) with degree $+1$ and $\lambda^{-1}$ multiplies every $E_{-\vec\alpha_j}$ with degree $-1$, so
     $$\mathcal A_\pm[\vec\phi_\vartheta](x;\lambda) = e^{\mp\vartheta}\mathcal A_\pm[\vec\phi](x'; e^\vartheta\lambda),$$
     and a boost is compensated by $\lambda\mapsto e^\vartheta\lambda$ alone.
   - In a general gradation with labels $(s_0,s_1,\dots,s_r)\ne(1,1,\dots,1)$ (such as the homogeneous gradation $(1,0,\dots,0)$), the spectral parameter $\tilde\lambda$ enters only on generators with $s_j > 0$, as $\tilde\lambda^{s_j}E_{\vec\alpha_j}$ and $\tilde\lambda^{-s_j}E_{-\vec\alpha_j}$, with total degree $N = \sum_{j=0}^r n_j s_j$. Under a boost, every term in $\mathcal A_\pm$ must still pick up the uniform factor $e^{\mp\vartheta}$. Rescaling $\tilde\lambda \mapsto e^{N\vartheta}\tilde\lambda$ multiplies $E_{\pm\vec\alpha_j}$ by $e^{\pm N s_j\vartheta}$ instead of $e^{\pm\vartheta}$. To supply the missing factor $e^{\pm(1 - N s_j)\vartheta}$, one must apply the constant Cartan gauge transformation $G(\vartheta) = \exp(\vartheta\vec\rho_s\cdot\vec H)$ where $\vec\alpha_j\cdot\vec\rho_s = 1 - N s_j$ for $j=1,\dots,r$ (and consequently $\vec\alpha_0\cdot\vec\rho_s = -\frac{1}{n_0}\sum_{j=1}^r n_j(1 - N s_j) = 1 - N s_0$), which is the grading difference between the principal and $(s_0,\dots,s_r)$ gradations. Thus `09-lagrangian.qmd:77` is **correct**.

3. **Spins of conserved charges (`09-lagrangian.qmd:84–88`) and explicit centralizer checks:**
   - Verified in Phase 1 Rank 10 (`check_rank10_exponents.py`) and Phase 2 (`check_phase2_items.py` Item 5):
     - $a_2^{(2)}$ ($k=2, h=3, kh=6$): exponents $1, 5\pmod 6$;
     - $a_4^{(2)}$ ($k=2, h=5, kh=10$): exponents $1, 3, 7, 9\pmod{10}$ (all odd integers not divisible by $2n+1=5$);
     - Dual pair $(c_2^{(1)}, d_3^{(2)})$: both have all odd integers $1, 3\pmod 4$;
     - Dual pair $(g_2^{(1)}, d_4^{(3)})$: $g_2^{(1)}$ ($k=1, h=6$) has $1, 5\pmod 6$, and $d_4^{(3)}$ ($k=3, h=4, kh=12$) has $1, 5, 7, 11\pmod{12}$, which are the exact same set of odd integers not divisible by $3$;
     - Dual pair $(f_4^{(1)}, e_6^{(2)})$: $f_4^{(1)}$ ($k=1, h=12$) has $1, 5, 7, 11\pmod{12}$ (all integers coprime to $6$), and $e_6^{(2)}$ ($k=2, h=9, kh=18$) has $1, 5, 7, 11, 13, 17\pmod{18}$ (also all integers coprime to $6$);
     - $a_n^{(1)}$ ($n\ge 2, k=1, h=n+1$): every integer not divisible by $h=n+1$.

4. **Spin-3 current (`eq-spin3-current`, `09-lagrangian.qmd:90–96`) and `@exr-spin3` (`09-lagrangian.qmd:168–170`):**
   - As found in Phase 1 Rank 8 (`check_ch09_lagrangian.py` Check 1):
     1. `eq-spin3-current` ($\partial_-\big[(\partial_+^2\phi)^2 + \frac{\beta^2}{4}(\partial_+\phi)^4\big] + \partial_+\big[m^2\cosh(\beta\phi)(\partial_+\phi)^2\big] = 0$) is **correct** for sinh-Gordon.
     2. However, at `09-lagrangian.qmd:96`, the sentence *"For sine-Gordon, replace $\beta^2$ by $-\beta^2$"* is **wrong**: replacing only $\beta^2\to -\beta^2$ while leaving $\cosh(\beta\phi)$ gives a non-zero divergence on the sine-Gordon equation $\partial_+\partial_-\phi = -\frac{m^2}{\beta}\sin(\beta\phi)$. One must replace $\beta\to i\beta$, which sends $\beta^2\to -\beta^2$ **and** $\cosh(\beta\phi)\to\cos(\beta\phi)$.
     3. Furthermore, with $\partial_\pm = \partial_t \pm \partial_x$, on a right-moving state of rapidity $\theta$ ($x\cosh\theta - t\sinh\theta$), $\partial_+$ multiplies by $m(\cosh\theta - \sinh\theta) = me^{-\theta}$, so the charge $\int dx\,(T_4 + \mathcal X_2)$ built from `eq-spin3-current` has Lorentz spin $-3$ and scales as $e^{-3\theta}$ (evaluating on the sine-Gordon kink to $-\frac{16m^3}{3\beta^2}e^{-3\theta}$), while exchanging $\partial_+\leftrightarrow\partial_-$ gives the spin-$+3$ charge scaling as $e^{+3\theta}$.

5. **Literature check of `[@feigin1993; @dgz1992; @niedermaier1994]` (`09-lagrangian.qmd:32`):**
   We checked all three primary sources in `/home/gustav/referee-copy/papers/`:
   - `feigin1993` (`hep-th_9310022.txt`, Feigin and Frenkel, *Integrable Systems and Quantum Groups*, 1993/1995): proves the existence of quantum higher-spin local integrals of motion for affine Toda theories by cohomological perturbation of the screening charges for general $\hat g$.
   - `dgz1992` (`hep-th_9201067.txt`, Delius, Grisaru and Zanon, *Nucl. Phys.* B382 (1992) 365): explicitly constructs renormalized quantum higher-spin conserved currents at one loop / all orders for non-simply-laced affine Toda theories (including $c_2^{(1)}, b_n^{(1)}, c_n^{(1)}$) and computes their anomalous threshold corrections.
   - `niedermaier1994` (`hep-th_9401078.txt`, Niedermaier, *Nucl. Phys.* B424 (1994) 184, "The quantum spectrum of conserved charges in affine Toda theories"): proves the existence and computes the eigenvalues of the full tower of quantum conserved charges in simply-laced affine Toda theories via the form-factor bootstrap.
   - Thus the citation `[@feigin1993; @dgz1992; @niedermaier1994]` at `09-lagrangian.qmd:32` is completely accurate.

### Verdict

- **Verdict:** **wrong** (`09-lagrangian.qmd:96` is wrong on two points: replacing only $\beta^2\to -\beta^2$ for sine-Gordon without $\cosh(\beta\phi)\to\cos(\beta\phi)$, and stating that `eq-spin3-current` built from $\partial_+ = \partial_t+\partial_x$ scales as $e^{+3\theta}$ rather than $e^{-3\theta}$; all other claims in `09-lagrangian.qmd:32–99` are **correct**).

### Paste-ready replacement text

- In `part3-classical/09-lagrangian.qmd`, line 96, replace:
  ```markdown
  The charge $Q_3=\int dx\,\big(T_4+\mathcal X_2\big)$ built from this current has spin 3, and on a one-particle state it takes the value $\chi\,m^3e^{3\theta}$, as in (-@eq-charge-eigenvalue). For sine-Gordon, replace $\beta^2$ by $-\beta^2$.
  ```
  with:
  ```markdown
  Because $\partial_+=\partial_t+\partial_x$ multiplies a mode of rapidity $\theta$ by $e^{-\theta}$, the charge $\int dx\,\big(T_4+\mathcal X_2\big)$ built from (-@eq-spin3-current) has spin $-3$ and scales as $e^{-3\theta}$; exchanging $\partial_+\leftrightarrow\partial_-$ gives the spin-$3$ charge $Q_3$, which on a one-particle state takes the value $\chi\,m^3e^{3\theta}$, as in (-@eq-charge-eigenvalue). For sine-Gordon, replace $\beta$ by $i\beta$, so $\beta^2\to-\beta^2$ and $\cosh(\beta\phi)\to\cos(\beta\phi)$.
  ```

---

## Item 6. Masses, three-point couplings, the charge basis and Dorey's rule (`09-lagrangian.qmd:102–161`)

### Overlap

- **Overlap with Phase 1:** **in both**.
  - **Covered in Phase 1:** Rank 11 (`09-lagrangian.qmd:122–138`, `eq-area-rule`, charge basis vs real orthonormal basis), Rank 12 (`09-lagrangian.qmd:142–160`, `prp-dorey-rule`, `fig-dorey-d4`, `exr-dorey-a2`), and Rank 14 (`09-lagrangian.qmd:114–117`, mass table for $a_n^{(1)}, d_n^{(1)}, e_8^{(1)}$ and sum rule).
  - **New in Phase 2:** Explicit numerical verification of the area rule (`eq-area-rule`) and Dorey's rule (`prp-dorey-rule`) — and their failure in a real orthonormal basis — for **$a_4^{(1)}$ and $e_6^{(1)}$** (`check_phase2_items.py` Item 6), supplementing $a_2^{(1)}$ and $d_4^{(1)}$ checked in Phase 1 (`check_ch09_lagrangian.py` Check 4).

### Check

1. **Masses table and sum rule (`09-lagrangian.qmd:104–120`):**
   - Verified in Phase 1 Rank 14 (`check_ch09_lagrangian.py` Check 3): the eigenvalues of $\mathsf M^2 = m^2\sum_{j=0}^r n_j\vec\alpha_j\vec\alpha_j^{\mathsf T}$ match the closed-form formulas in the table (`09-lagrangian.qmd:114–116`) for $a_n^{(1)}$, $d_n^{(1)}$ ($n=4,5,6$), and all 8 trigonometric expressions for $e_8^{(1)}$ to $8.88\times 10^{-16}$, and satisfy $\sum_{a=1}^r m_a^2 = m^2\sum_{j=0}^r n_j\vec\alpha_j^2 = 2hm^2$.

2. **Area rule (`eq-area-rule`), Dorey's rule (`prp-dorey-rule`), and the charge eigenbasis vs real orthonormal basis for $a_2^{(1)}, a_4^{(1)}, d_4^{(1)}, e_6^{(1)}$ (`09-lagrangian.qmd:130–160`):**
   - In `check_ch09_lagrangian.py` (Check 4) and `check_phase2_items.py` (Item 6), we computed $C_{abc} = m^2\beta\sum_{j=0}^r n_j(\vec\alpha_j\cdot\vec e_a)(\vec\alpha_j\cdot\vec e_b)(\vec\alpha_j\cdot\vec e_c)$ ($m=\beta=1$) in the centre-character eigenbasis (normalized by $\vec e_a\cdot\vec e_{\bar a} = 1$) and compared $\lvert C_{abc}\rvert$ with $\frac{4\beta}{\sqrt{h}}\Delta_{abc}$:
     ```text
     a_2^(1) (h=3) in Z_3 charge basis: |C_111|=|C_222|=3.000000, (4/sqrt(3))Delta_111=3.000000, diff=4.44e-16; C_112=C_122=0
     a_4^(1) (h=5) in Z_5 charge basis: max ||C_abc| - (4/sqrt(5)) Delta_abc| = 8.24e-16
     d_4^(1) (h=6) in Z_2 x Z_2 charge basis: max ||C_abc| - (4/sqrt(6)) Delta_abc| = 2.22e-15
     e_6^(1) (h=12) in Z_3 charge basis: 18 non-zero couplings, max ||C_abc| - (4/sqrt(12)) Delta_abc| = 3.55e-15
     ```
   - **Check in a real orthonormal basis for complex-conjugate pairs ($a \ne \bar a$):**
     - In $a_2^{(1)}$: in the real orthonormal basis $\vec u_1 = (\vec e_1+\vec e_2)/\sqrt{2}$, $\vec u_2 = (\vec e_1-\vec e_2)/(i\sqrt{2})$, we find $\lvert C_{u_1 u_1 u_1}\rvert = \frac{3}{\sqrt{2}} \ne 3 = \frac{4}{\sqrt{3}}\Delta_{111}$, violating `eq-area-rule`.
     - In $a_4^{(1)}$: in the real orthonormal basis $\vec u_1 = (\vec e_1+\vec e_4)/\sqrt{2}$, $\vec u_2 = (\vec e_2+\vec e_3)/\sqrt{2}$, we find $\lvert C_{u_1 u_1 u_1}\rvert = 0$ and $\lvert C_{u_1 u_1 u_2}\rvert = 0 \ne 1.175571 = \frac{4}{\sqrt{5}}\Delta_{112}$, while $\lvert C_{u_1 v_1 u_2}\rvert$ and $\lvert C_{v_1 v_1 v_2}\rvert$ differ from $\frac{4}{\sqrt{5}}\Delta_{112}$.
     - In $e_6^{(1)}$: for the conjugate pair with $m_a = 2.175328$, in the real orthonormal basis $\vec u_a = (\vec e_a+\vec e_{\bar a})/\sqrt{2}$ we find $\lvert C_{u_a u_a u_a}\rvert = 0.978591 \ne 2.366025 = \frac{4}{\sqrt{12}}\Delta_{aaa}$.
   - **What is wrong in `09-lagrangian.qmd:130` for even $d_{2m}^{(1)}$ (found in Phase 1 Rank 11):**
     Line 130 states:
     > *"For degenerate masses ($a_n$, $d_n$, $e_6$) the $\vec e_a$ must be charge eigenvectors, the complex eigenvectors of the centre of the affine diagram (of $\mathbb Z_h$ for $a_n$), normalized by $\vec e_a\cdot\vec e_{\bar a}=1$ so that $\phi_{\bar a}=\bar\phi_a$. Particle $a$ is the eigenvector with the centre character of node $a$. In a real orthonormal basis of a degenerate pair neither (-@eq-area-rule) nor @prp-dorey-rule holds."*
     For $d_{2m}^{(1)}$ with **even** $n=2m$ (such as $d_4^{(1)}$, illustrated in `@fig-dorey-d4` on the very same page!), the centre is $\mathbb Z_2\times\mathbb Z_2$, all of whose characters are **real** ($\pm 1$), every particle is self-conjugate ($\bar a = a$), and the charge eigenbasis of the degenerate spinor pair $(s, s')$ (and of the outer triality triplet $(1, s, s')$ in $d_4^{(1)}$) **is** a real orthonormal basis ($\vec e_a\in\mathbb R^r$, $\vec e_a\cdot\vec e_b = \delta_{ab}$) in which both `eq-area-rule` and `@prp-dorey-rule` hold! It is only for complex-conjugate degenerate pairs ($a\ne\bar a$, namely $a_n^{(1)}$ with $a\ne h/2$, $d_{2m+1}^{(1)}$ spinors, and $e_6^{(1)}$) that the charge eigenvectors are complex and any real orthonormal basis violates `eq-area-rule` and `@prp-dorey-rule`.

3. **Dorey's rule (`09-lagrangian.qmd:142–158`, `prp-dorey-rule`, `fig-dorey-d4`):**
   - Verified in Phase 1 Rank 12 (`check_ch09_lagrangian.py` Check 4):
     - If $a$ is white, $s_{\rm black}\vec\lambda_a = \vec\lambda_a$ and $s_{\rm white}\vec\lambda_a = \vec\lambda_a - \vec\alpha_a$, so $\vec\gamma_a = (1 - s_{\rm white}s_{\rm black})\vec\lambda_a = \vec\alpha_a$.
     - If $a$ is black, $s_{\rm white}\vec\lambda_a = \vec\lambda_a$ and $s_{\rm black}\vec\lambda_a = \vec\lambda_a - \vec\alpha_a$, so $\vec\gamma_a = (1 - s_{\rm white}s_{\rm black})\vec\lambda_a = \vec\lambda_a - s_{\rm white}(\vec\lambda_a - \vec\alpha_a) = s_{\rm white}\vec\alpha_a$.
     - In $d_4^{(1)}$ ($h=6$), the four Coxeter orbits $\Omega_1, \Omega_c, \Omega_s, \Omega_{s'}$ have 6 roots each, with projected radii in the Coxeter plane in the ratio $1 : \sqrt{3} : 1 : 1$, and two roots from $\Omega_1$ and one root from $\Omega_c$ close a triangle $\vec\gamma_1 + \vec\gamma_1' + \vec\gamma_c = 0$ with side ratios $1, 1, \sqrt{3}$ and interior angles $\pi/6, \pi/6, 2\pi/3$, matching `@fig-dorey-d4`.

### Verdict

- **Verdict:** **wrong** (`09-lagrangian.qmd:130` overclaims that the charge eigenvectors of $d_n$ are always complex and that `eq-area-rule` and `@prp-dorey-rule` fail in any real orthonormal basis of a degenerate pair, whereas for even $d_{2m}^{(1)}$ the centre characters are $\pm 1$ and the charge eigenbasis of $(s,s')$ is itself real orthonormal; all other claims in `09-lagrangian.qmd:102–161` are **correct**).

### Paste-ready replacement text

- In `part3-classical/09-lagrangian.qmd`, line 130, replace:
  ```markdown
  For degenerate masses ($a_n$, $d_n$, $e_6$) the $\vec e_a$ must be charge eigenvectors, the complex eigenvectors of the centre of the affine diagram (of $\mathbb Z_h$ for $a_n$), normalized by $\vec e_a\cdot\vec e_{\bar a}=1$ so that $\phi_{\bar a}=\bar\phi_a$. Particle $a$ is the eigenvector with the centre character of node $a$. In a real orthonormal basis of a degenerate pair neither (-@eq-area-rule) nor @prp-dorey-rule holds.
  ```
  with:
  ```markdown
  For degenerate masses ($a_n$, $d_n$, $e_6$) the $\vec e_a$ must be charge eigenvectors, the eigenvectors of the centre of the affine diagram (of $\mathbb Z_h$ for $a_n$; see @sec-hirota-form-dn1-en1), normalized by $\vec e_a\cdot\vec e_{\bar a}=1$ so that $\phi_{\bar a}=\bar\phi_a$ (which are complex for conjugate pairs $a\ne\bar a$ in $a_n$, $d_{2m+1}$ and $e_6$, and real orthonormal for even $d_{2m}$, where $\bar a=a$). Particle $a$ is the eigenvector with the centre character of node $a$. In a real orthonormal basis of a complex-conjugate pair ($a\ne\bar a$), neither (-@eq-area-rule) nor @prp-dorey-rule holds.
  ```

---

## Item 7. Solitons, topological charges, Lie duality and Part VIII pointers (`10-solitons.qmd`)

### Overlap

- **Overlap with Phase 1:** **in both**.
  - **Covered in Phase 1:** Rank 5 (`10-solitons.qmd:22–25, 104, 158–160`, Lie duality and `kneipp1996`), Rank 6 (`10-solitons.qmd:76–77`, polynomial degrees of folded solitons and $c_2^{(1)}$ static $\tau_j$), Rank 9 (`10-solitons.qmd:82–92, 154–156`, topological charges of $a_n^{(1)}$ and `mcghee1994`), Rank 15 (`10-solitons.qmd:10`, embedded sine-Gordon kinks), Rank 17 (`10-solitons.qmd:68–79`, vertex operators), and Rank 19 (`10-solitons.qmd:126, 162–164`, fusion limit).
  - **New in Phase 2:**
    1. Explicit check of folded polynomial degrees for $g_2^{(1)}$ folded from $d_4^{(1)}$ (`check_phase2_items.py` Item 7).
    2. Explicit verification of the Part VIII cross-references `@prp-stored-charge` (`39-defect-solitons.qmd:76`) and `@prp-bianchi` (`37-defect-sine-gordon.qmd:43`).

### Check

1. **Hirota bilinear equations (`eq-hirota-bilinear`) and $a_n^{(1)}$ single soliton (`eq-an-soliton`, `10-solitons.qmd:31–64`):**
   - Substituting $\vec\phi = \frac{i}{\beta}\sum_j \frac{2\vec\alpha_j}{\alpha_j^2}\ln\tau_j$ into `eq-toda-eom` ($\partial_+\partial_-\vec\phi = \frac{im^2}{\beta}\sum_k n_k\vec\alpha_k e^{i\beta\vec\alpha_k\cdot\vec\phi}$) and using $e^{i\beta\vec\alpha_k\cdot\vec\phi} = \prod_j \tau_j^{-a_{jk}}$ and $\sum_j n_j\vec\alpha_j = 0$ gives $\partial_+\partial_-\ln\tau_j = m^2 n_j\frac{\vec\alpha_j^2}{2}\big(\prod_l \tau_l^{-a_{lj}} - 1\big)$. Since $\partial_+\partial_-\ln\tau_j = (\partial_t^2 - \partial_x^2)\ln\tau_j = -\frac12(D_x^2 - D_t^2)\tau_j\cdot\tau_j / \tau_j^2$, this gives `eq-hirota-bilinear`:
     $$\tfrac12(D_x^2 - D_t^2)\tau_j\cdot\tau_j = m^2 n_j\frac{\vec\alpha_j^2}{2}\Big(\tau_j^2 - \prod_{l\ne j}\tau_l^{-a_{lj}}\Big).$$
   - For $a_n^{(1)}$ (`eq-an-soliton`, $\tau_j = 1 + \omega^{ja}E$), $\frac12(D_x^2 - D_t^2)\tau_j\cdot\tau_j = m_a^2(\cosh^2\theta - \sinh^2\theta)\omega^{ja}E = m_a^2\omega^{ja}E$, while $\tau_j^2 - \tau_{j-1}\tau_{j+1} = (2 - \omega^a - \omega^{-a})\omega^{ja}E = 4\sin^2\frac{\pi a}{h}\,\omega^{ja}E$, which matches for $m_a = 2m\sin\frac{\pi a}{h}$.

2. **Embedded sine-Gordon kinks (`10-solitons.qmd:10`):**
   - Verified in Phase 1 Rank 15 (`check_ch10_solitons.py` Check 4): for species $a=n$ of $a_{2n-1}^{(1)}$ ($h=2n$) at $\operatorname{Im}\xi = \pi/2$, $\omega^{jn} = (-1)^j$ and $E = iR$ ($R > 0$), so $\tau_{2k} = 1 + iR$ and $\tau_{2k+1} = 1 - iR = \bar\tau_{2k}$, giving the strictly real field $\vec\phi(x,t) = -\frac{2}{\beta}\big(\sum_{k=0}^{n-1}\vec\alpha_{2k}\big)\arctan R$.

3. **Polynomial degrees of folded solitons and the $c_2^{(1)}$ and $g_2^{(1)}$ examples (`10-solitons.qmd:76–77`):**
   - **Check of $c_2^{(1)}$ (`check_ch10_solitons.py` Check 2):** Substituting $\tau_0 = 1 + 2E + \frac12 E^2$, $\tau_1 = 1 + \frac12 E^2$, $\tau_2 = 1 - 2E + \frac12 E^2$ with $E = \exp(\sqrt{2}mx)$ into `eq-hirota-bilinear` for $c_2^{(1)}$ ($n_0=1, n_1=2, n_2=1$, $\vec\alpha_0^2=\vec\alpha_2^2=2, \vec\alpha_1^2=1$, so $n_j\vec\alpha_j^2/2 = 1$ for all $j$, with $\tau_0^2 - \tau_1^2 = 4E(1+\frac12 E^2)$, $\tau_1^2 - \tau_0\tau_2 = 4E^2$, $\tau_2^2 - \tau_1^2 = -4E(1+\frac12 E^2)$) gives zero residual at every node $j=0,1,2$.
   - **Check of $g_2^{(1)}$ folded from $d_4^{(1)}$ (`check_phase2_items.py` Item 7):** In $d_4^{(1)}$, the outer nodes $0,1,3,4$ have $n_j^{\rm parent} = 1$ and the central node $2$ has $n_2^{\rm parent} = 2$. Folding $\{1,3,4\}$ to the short node of $g_2^{(1)}$ gives $n_j = (1,2,3)$ and $n_j^\vee = n_j\vec\alpha_j^2/2 = (1,2,1) = n_j^{\rm parent}$. The single central soliton ($\lvert O\rvert = 1$) has degrees $\lvert O\rvert n_j^\vee = (1,2,1)$, and the coincident triple of outer solitons ($\lvert O\rvert = 3$) has degrees $\lvert O\rvert n_j^\vee = (3,6,3)$.
   - **Exception at $a_{2n}^{(2)}$ and notation slip $\tau_1 = \tau_3$ (found in Phase 1 Rank 6):**
     As shown in Phase 1 Rank 6, for all 8 non-simply-laced foldings in `@tbl-foldings`, $n_j^\vee = k n_j\vec\alpha_j^2/2$ equals the parent Kac label $n_j^{\rm parent}$ (`@hyp-D` in `18-jordan-chains.qmd:323`). However, for the self-dual twisted algebra $a_{2n}^{(2)}$ folded from $a_{2n}^{(1)}$ (`06-affine-algebras.qmd:113`), the folded node $\{n,n+1\}$ joins adjacent roots so $\vec\alpha_n^2 = 1/2$, $n_j^\vee = (2,\dots,2,1) \ne n_j^{\rm parent} = (1,\dots,1)$, and a coincident pair ($\lvert O\rvert = 2$) has degree $\lvert O\rvert n_j^{\rm parent} = 2$ at every node (`18-jordan-chains.qmd:350`), not $\lvert O\rvert n_j^\vee = 4$. Also, in $c_2^{(1)}$ there are only nodes $0,1,2$; writing $\tau_1 = \tau_3 = 1 + \frac12 E^2$ uses the parent $a_3^{(1)}$ node index $3$ without explanation.

4. **Topological charges and `[@mcghee1994]` (`10-solitons.qmd:80–93`):**
   - **Exact source check (`hep-th_9307035.txt`, McGhee, *Int. J. Mod. Phys.* A9 (1994) 2645):**
     McGhee explicitly proves the $h/\gcd(a,h)$ formula for all $a_n^{(1)}$! In the Abstract and Introduction (`hep-th_9307035.txt`, lines 18–22 and 82–88):
     > *"First the $a_n^{(1)}$ series is investigated and we find that for the $a^{\rm th}$ soliton ($a=1,\dots,n$) on the $a^{\rm th}$ node of the Dynkin diagram, the number of topological charges obtained is fewer than dimensionality of the $a^{\rm th}$ fundamental representation whenever $a$ and $n+1$ are not coprime. In fact the number of distinct charges found is $\tilde h_a = \frac{h}{\gcd(a,h)}$ where $h=n+1$ is the Coxeter number and $\gcd(a,h)$ is the greatest common divisor of $a$ and $h$."*
     And in Section 3 (`hep-th_9307035.txt`, lines 435–450), McGhee proves the jumping rule as $\operatorname{Im}\xi$ varies by $2\pi$.
   - **Sign and representation check (`check_ch10_solitons.py` Check 1):** Verified in Phase 1 Rank 9 that $\frac{\beta}{2\pi}\vec\phi(+\infty) = -\frac{1}{2\pi}\sum_{j=0}^{h-1}\vec\alpha_j[\vartheta_j]_{(-\pi,\pi)}$ is a weight of $\Lambda^a\mathbb C^h$ (not $\Lambda^{h-a}\mathbb C^h$), jumping at $h/\gcd(a,h)$ values of $\operatorname{Im}\xi \in [0,2\pi)$, $\gcd(a,h)$ phases at once.

5. **Soliton masses (`eq-soliton-mass`), `@exr-soliton-mass` hint (`10-solitons.qmd:150–152`), and Lie duality / `[@kneipp1996]` (`10-solitons.qmd:104`):**
   - **Check of `@exr-soliton-mass` (`10-solitons.qmd:150–152`):**
     For a static solution $\vec\phi(x)$, the equation of motion (`eq-toda-eom` with $\partial_+\partial_- = -\partial_x^2$) is $\vec\phi'' = \nabla V_{\rm imag}(\vec\phi)$. Multiplying by $\vec\phi'$ (without complex conjugation) and integrating gives the first integral $\frac12\vec\phi'\cdot\vec\phi' = V_{\rm imag}(\vec\phi)$ (since $\vec\phi'\to 0$ and $V_{\rm imag}\to 0$ as $x\to -\infty$). Thus the static energy is
     $$E = \int_{-\infty}^\infty \Big(\tfrac12\vec\phi'\cdot\vec\phi' + V_{\rm imag}(\vec\phi)\Big)\,dx = 2\int_{-\infty}^\infty V_{\rm imag}(\vec\phi)\,dx.$$
     From the intermediate Hirota equation (`10-solitons.qmd:42`, $-\partial_x^2\ln\tau_j = m^2 n_j\frac{\vec\alpha_j^2}{2}\big(\prod_l \tau_l^{-a_{lj}} - 1\big) = m^2 n_j\frac{\vec\alpha_j^2}{2}\big(e^{i\beta\vec\alpha_j\cdot\vec\phi} - 1\big)$), dividing by $\beta^2\vec\alpha_j^2/2$ and summing over $j$ gives
     $$V_{\rm imag}(\vec\phi) = -\frac{m^2}{\beta^2}\sum_j n_j\big(e^{i\beta\vec\alpha_j\cdot\vec\phi} - 1\big) = \frac{1}{\beta^2}\sum_j \frac{2}{\vec\alpha_j^2}\partial_x^2\ln\tau_j.$$
     Integrating $2V_{\rm imag}$ over $x\in(-\infty,\infty)$ gives $E = \frac{2}{\beta^2}\sum_j \frac{2}{\vec\alpha_j^2}\big[\partial_x\ln\tau_j\big]_{-\infty}^\infty$. For a simply-laced static soliton of species $a$ ($\vec\alpha_j^2 = 2$, $\tau_j$ of degree $n_j$ in $E = e^{m_a x + \xi}$), $\partial_x\ln\tau_j \to n_j m_a$ as $x\to +\infty$ and $0$ as $x\to -\infty$, so $E = \frac{2}{\beta^2}\big(\sum_{j=0}^r n_j\big)m_a = \frac{2hm_a}{\beta^2} = M_a$. Every step of the hint in `@exr-soliton-mass` is **correct**.
   - **Exact source check of `[@kneipp1996]` (`hep-th_9404030.txt`, Kneipp and Olive, *Nucl. Phys.* B408 (1993) 565 — note that in `references.bib` the key `kneipp1996` is actually a 1995/1996 paper or `kneipp1993`? Wait: let's check `references.bib`!):**
     In `hep-th_9404030.txt` (Kneipp and Olive, Section 6, lines 796–801 and 945–986), Kneipp and Olive treat **both** twisted foldings ($\tau\in W_0(g)$, $M_{\langle i\rangle} = M_i$, eq. (6.1)) **and** untwisted non-simply-laced foldings ($\tau\in\operatorname{Aut}\Delta(g)$, $M_{\langle i\rangle} = \lvert\langle i\rangle\rvert M_i = \frac{4h}{\alpha_{\langle i\rangle}^2\beta^2}m_{\langle i\rangle}$, eq. (6.9) and lines 951–986), explicitly deriving the soliton masses of untwisted non-simply-laced theories by folding the energy–momentum tensor and confirming Olive–Turok–Underwood (`hep-th/9305160`, Section 7)!

6. **Part VIII pointers `@prp-stored-charge` (`10-solitons.qmd:90`) and `@prp-bianchi` (`10-solitons.qmd:120`):**
   - `@prp-stored-charge` (`part8-defects/39-defect-solitons.qmd:76–79`): proves that for $a_n^{(1)}$ a species-$a$ soliton transmitted through a type-I defect experiences a phase shift $\tau_j \mapsto \tau_j$ with $\omega^{ja}E \mapsto d_{j,a}\omega^{ja}E$, which shifts $\operatorname{Im}\xi$ and moves the soliton's topological charge along the same cycle of $h/\gcd(a,h)$ weights without ever producing a missing charge. Matches `10-solitons.qmd:90`.
   - `@prp-bianchi` (`part8-defects/37-defect-sine-gordon.qmd:43–49`): proves Bianchi's permutability theorem $\tan\frac{\beta_{\rm SG}(\phi_{12}-\phi_0)}{4} = \frac{\kappa_1+\kappa_2}{\kappa_1-\kappa_2}\tan\frac{\beta_{\rm SG}(\phi_1-\phi_2)}{4}$ with $\kappa_i = e^{-\theta_i}$, whose square of the inverse coefficient is $\big(\frac{\kappa_1-\kappa_2}{\kappa_1+\kappa_2}\big)^2 = \tanh^2\frac{\theta_1-\theta_2}{2} = A_{11}(\theta_1-\theta_2)$. Matches `10-solitons.qmd:120`.

### Verdict

- **Verdict:** **wrong** on two items already found in Phase 1:
  1. `10-solitons.qmd:76–77` (Phase 1 Rank 6): states that a folded soliton of orbit $O$ has degree $\lvert O\rvert n_j^\vee$ for all folded theories (which fails for $a_{2n}^{(2)}$ where the degree is $\lvert O\rvert n_j^{\rm parent} = 2$ while $\lvert O\rvert n_j^\vee = 4$), and writes $\tau_1 = \tau_3$ in $c_2^{(1)}$ using the parent $a_3^{(1)}$ node index $3$.
  2. `10-solitons.qmd:126` and `10-solitons.qmd:162–164` (`#exr-fusion`, Phase 1 Rank 19): omits the condition $a+b < h$ for $m_a e^{i\pi b/h} + m_b e^{-i\pi a/h} = m_{a+b}$ (for $a+b > h$, $\sin\frac{\pi(a+b)}{h} < 0$ and the fused species $a+b-h$ arises at $\theta_1 = \theta - i\pi(h-b)/h, \theta_2 = \theta + i\pi(h-a)/h$) and states that $\theta_1-\theta_2 = i\pi(a+b)/h$ is fixed before sending $\xi_1,\xi_2\to -\infty$ rather than as a joint limit.
  3. `10-solitons.qmd:128–136` (`#sec-breathers`, Phase 1 Rank 3): see Item 8 below.

### Paste-ready replacement text

- See Phase 1 Rank 6 above for `10-solitons.qmd:76–77` and Phase 1 Rank 19 above for `10-solitons.qmd:126` and `10-solitons.qmd:162–164`:
  - In `part3-classical/10-solitons.qmd`, lines 76–77, replace:
    ```markdown
    For folded theories the relevant label is $n_j^\vee$, the level of $L(\Lambda_j)$. A coincident orbit $O$ of $\lvert O\rvert$ parent solitons has degree $\lvert O\rvert n_j^\vee$. In $c_n^{(1)}$ every $n_j^\vee$ is 1, so a coincident pair of species $a$ and $2n-a$ has degree 2 at every node, including the long end nodes, whose Kac label is 1, and the species-$n$ soliton has degree 1 at every node, including the short nodes, whose Kac label is 2. For $c_2^{(1)}$, folded from $a_3^{(1)}$, the static pair of species 1 and 3 has $\tau_0=1+2E+\frac12E^2$, $\tau_1=\tau_3=1+\frac12E^2$ and $\tau_2=1-2E+\frac12E^2$.
    ```
    with:
    ```markdown
    For theories folded from a simply-laced parent (@tbl-foldings), the relevant label at node $j$ is the parent Kac label $n_j^{\rm parent}$, which for every untwisted and twisted pair in @tbl-foldings equals $n_j^\vee$, the level of $L(\Lambda_j)$ (for $a_{2n}^{(2)}$ folded from $a_{2n}^{(1)}$, every $n_j^{\rm parent}=1$). A coincident orbit $O$ of $\lvert O\rvert$ parent solitons has degree $\lvert O\rvert n_j^{\rm parent}$. In $c_n^{(1)}$ every $n_j^\vee=n_j^{\rm parent}$ is 1, so a coincident pair of species $a$ and $2n-a$ has degree 2 at every node, including the long end nodes, whose Kac label is 1, and the species-$n$ soliton has degree 1 at every node, including the short nodes, whose Kac label is 2. For $c_2^{(1)}$, folded from $a_3^{(1)}$, the static pair of species 1 and 3 has $\tau_0=1+2E+\frac12E^2$, $\tau_1=1+\frac12E^2$ (from $\tau_1=\tau_3$ in $a_3^{(1)}$) and $\tau_2=1-2E+\frac12E^2$.
    ```
  - In `part3-classical/10-solitons.qmd`, line 126 and lines 162–164, replace:
    ```markdown
    **Fusion.** At $\theta_1-\theta_2=i\pi(a+b)/h$, the denominator of $A_{ab}$ vanishes. Take $\theta_1=\theta+i\pi b/h$ and $\theta_2=\theta-i\pi a/h$. Then $E_1E_2$ has the exponent of a single soliton of species $a+b$ with rapidity $\theta$, because $m_ae^{i\pi b/h}+m_be^{-i\pi a/h}=m_{a+b}$, and the phases combine to $\omega^{j(a+b)}$. Letting $\xi_1,\xi_2\to-\infty$ while keeping $A_{ab}e^{\xi_1+\xi_2}$ fixed reduces (-@eq-n-soliton) to the single soliton of species $a+b$.
    ```
    with:
    ```markdown
    **Fusion.** For $a+b<h$, the denominator of $A_{ab}$ vanishes at $\theta_1-\theta_2=i\pi(a+b)/h$. Taking $\theta_1\to\theta+i\pi b/h$ and $\theta_2\to\theta-i\pi a/h$, $E_1E_2$ acquires the exponent of a single soliton of species $a+b$ with rapidity $\theta$, because $m_ae^{i\pi b/h}+m_be^{-i\pi a/h}=m_{a+b}$, and the phases combine to $\omega^{j(a+b)}$. Sending $\xi_1,\xi_2\to-\infty$ jointly with $\theta_1-\theta_2\to i\pi(a+b)/h$ so that $A_{ab}(\theta_1-\theta_2)e^{\xi_1+\xi_2}\to e^\xi$ reduces (-@eq-n-soliton) to the single soliton of species $a+b$ (for $a+b>h$, take $\theta_1\to\theta-i\pi(h-b)/h$ and $\theta_2\to\theta+i\pi(h-a)/h$ to obtain species $a+b-h$).
    ```
    and in `@exr-fusion` (`10-solitons.qmd:162–164`), replace:
    ```markdown
    ::: {#exr-fusion}
    Verify $m_ae^{i\pi b/h}+m_be^{-i\pi a/h}=m_{a+b}$ for $m_a=\sin\frac{\pi a}h$, and carry out the fusion limit of @sec-multi-solitons for two solitons of $a_2^{(1)}$.
    :::
    ```
    with:
    ```markdown
    ::: {#exr-fusion}
    For $a+b<h$, verify $m_ae^{i\pi b/h}+m_be^{-i\pi a/h}=m_{a+b}$ for $m_a=\sin\frac{\pi a}h$, and carry out the joint fusion limit of @sec-multi-solitons for two species-1 solitons of $a_2^{(1)}$.
    :::
    ```

---

## Item 8. Singular solutions, breathers and the missing-charge problem (`10-solitons.qmd:128–139`, `11-reality.qmd:58–75`)

### Overlap

- **Overlap with Phase 1:** **in both**.
  - **Covered in Phase 1:** Rank 3 (`10-solitons.qmd:128–136`, `11-reality.qmd:28, 87`, breather parameter range $u$ and singularity condition), Rank 7 (`11-reality.qmd:58–62, 90–92`, zeros of $\tau_j$, double poles with zero residue, contour deformation vs jump of topological charge), and Rank 20 (`11-reality.qmd:62, 64–75`, `10-solitons.qmd:138, 181`, Takács–Watts `tw1999` and missing-charge cross-references).

### Check

1. **Zeros of $\tau_j$ in the real $(x,t)$ plane, double poles with zero residue, and whether integrated charges are unaffected (`11-reality.qmd:60–61`, `exr-singular-points`):**
   - Verified in Phase 1 Rank 7 (`check_ch11_reality_breathers.py` Check 5):
     1. *Double pole with zero residue:* Near a zero $x_0(t)$ of $\tau_j(\cdot,t)$ of multiplicity $m_0\ge 1$, $\tau_j(x,t) = (x-x_0)^{m_0}g(x)$ with $g(x_0)\ne 0$, so $\partial_x\ln\tau_j = \frac{m_0}{x-x_0} + \frac{g'}{g}$ and
        $$T_{00} = \frac{2}{\beta^2}\sum_j \frac{2}{\vec\alpha_j^2}\partial_x^2\ln\tau_j = \frac{2}{\beta^2}\sum_j \frac{2}{\vec\alpha_j^2}\left(-\frac{m_0}{(x-x_0)^2} + \partial_x\frac{g'}{g}\right).$$
        Because $\partial_x^2\ln\tau_j$ is an exact $x$-derivative of a function whose only pole is simple, its $(x-x_0)^{-1}$ residue is **identically zero**.
     2. *Three corrections required in `11-reality.qmd:60–61` and `exr-singular-points` (`11-reality.qmd:90–92`):*
        - (a) For a single soliton at *real* rapidity $\theta$, $\tau_j(x,t)$ depends on $(x,t)$ only through the single real linear combination $x\cosh\theta - t\sinh\theta$, so the two real equations $\operatorname{Re}\tau_j = \operatorname{Im}\tau_j = 0$ are functions of one real variable: its real zero set is **never** an isolated point in $\mathbb R^2$ (it is empty for generic $\operatorname{Im}\xi$ and a straight line for tuned $\operatorname{Im}\xi$). Zeros are generically isolated points in $\mathbb R^2$ only when $x$ and $t$ enter with $\mathbb R$-linearly independent coefficients (e.g. at complex rapidity $\operatorname{Im}\theta\notin\pi\mathbb Z$ or in multi-soliton collisions with distinct rapidities).
        - (b) A double pole $-m_0/(x-x_0)^2$ on the real $x$-axis is **not** integrable on $\mathbb R$ in the Lebesgue, improper Riemann, or Cauchy principal-value sense ($\int_{\mathbb R} \frac{dx}{(x-x_0)^2} = +\infty$). Local conserved charges are unaffected only when defined by **deforming the integration contour into the complex $x$-plane** around $x_0$ (which is independent of whether the contour passes above or below $x_0$ precisely because the residue vanishes).
        - (c) The **topological charge** $\vec Q_{\rm top} = -\frac{1}{2\pi}\sum_j \vec\alpha_j^\vee[\operatorname{Im}\ln\tau_j]_{-\infty}^{+\infty}$ defined along the real $x$-axis is **not** unaffected when an isolated zero $x_0(t)$ of $\tau_j$ crosses the real axis at $t=t_0$: the winding phase $[\operatorname{Im}\ln\tau_j]_{-\infty}^{+\infty}$ jumps by $\pm 2\pi$, so $\vec Q_{\rm top}$ jumps by the coroot $\mp\vec\alpha_j^\vee$!

2. **Breathers (`10-solitons.qmd:128–136`, `#fig-breather-theta` at `11-reality.qmd:28`, and `#exr-breather-real` at `11-reality.qmd:87`):**
   - Verified in Phase 1 Rank 3 (`check_ch11_reality_breathers.py` Check 3) against Harder, Iskandar and McGhee (`hep-th/9409035`, `harder1995`, Section 3, eqs. (3.5)–(3.16)):
     - For a soliton–antisoliton pair $(a, \bar a = h-a)$ of $a_n^{(1)}$ at $\theta_1 = \theta_0 + iu$, $\theta_2 = \theta_0 - iu$ ($u \in (0,\pi/2)$), the interaction coefficient is
       $$A_{a,h-a}(2iu) = \frac{\cos 2u - \cos\frac{2\pi a}{h}}{\cos 2u - \cos\pi} = 1 - \frac{\sin^2\frac{\pi a}{h}}{\cos^2 u}.$$
     - When $a = h/2$ (self-conjugate, as in sine-Gordon $a_1^{(1)}$ where $h=2, a=1$), $A_{a,a}(2iu) = -\tan^2 u < 0$ for **all** $u \in (0,\pi/2)$.
       - If $\xi_2 = \bar\xi_1$: then $E_2 = \bar E_1$, $E_1 E_2 = \lvert E_1\rvert^2 > 0$, so $\tau_j = 1 + 2\operatorname{Re}(\omega^{ja}E_1) - \tan^2 u\,\lvert E_1\rvert^2$ is **real** ($\bar\tau_j = \tau_j$) and changes sign from $+1$ as $x\to -\infty$ to $-\infty$ as $x\to +\infty$, vanishing along curves in the $(x,t)$ plane!
       - If $\xi_2 = \bar\xi_1 + i\pi$: then $E_2 = -\bar E_1$, $E_1 E_2 = -\lvert E_1\rvert^2 < 0$, so $A_{11}(2iu)E_1 E_2 = +\tan^2 u\,\lvert E_1\rvert^2 > 0$. For $a_1^{(1)}$ ($\omega = -1$), writing $E_1 = e^{X - iT + \xi_1}$ with $X = m_1\cos u\,(x\cosh\theta_0 - t\sinh\theta_0)$ and $T = m_1\sin u\,(t\cosh\theta_0 - x\sinh\theta_0)$, and choosing $\operatorname{Re}\xi_1 = -\ln\tan u$, $\operatorname{Im}\xi_1 = 0$ ($\xi_2 = i\pi$), gives
         $$\tau_0 = 1 + E_1 - \bar E_1 + (\tan^2 u)\lvert E_1\rvert^2 = 1 + e^{2X} - 2i\cot u\,e^X\sin T,$$
         $$\tau_1 = 1 - E_1 + \bar E_1 + (\tan^2 u)\lvert E_1\rvert^2 = 1 + e^{2X} + 2i\cot u\,e^X\sin T = \bar\tau_0,$$
         so $\overline{\tau_0/\tau_1} = \tau_1/\tau_0$, $\operatorname{Re}\tau_j = 1 + e^{2X} \ge 1 > 0$ (nowhere zero!), and with $u_{\rm SG} = \frac\pi2 - u$, $m_0 = m_1 = 2m$, $\beta_{\rm SG} = \sqrt{2}\beta$:
         $$\phi = \frac{i\sqrt{2}}{\beta}\ln\frac{\tau_0}{\tau_1} = \frac{2i}{\beta_{\rm SG}}\ln\frac{\tau_0}{\bar\tau_0} = \frac{4}{\beta_{\rm SG}}\arctan\left(\frac{\tan u_{\rm SG}\sin(m_0 t\cos u_{\rm SG})}{\cosh(m_0 x\sin u_{\rm SG})}\right),$$
         recovering `eq-sg-breather` (`03-sine-gordon.qmd:69`)!
     - **What goes wrong for $a \ne h/2$ ($h \ge 3$, e.g. $a_2^{(1)}$):**
       For $1 \le a < h/2$, $A_{a,h-a}(2iu) = 1 - \frac{\sin^2(\pi a/h)}{\cos^2 u}$ **changes sign at $u = \frac\pi2 - \frac{\pi a}{h} \in (0,\frac\pi2)$**:
       - For $\frac\pi2 - \frac{\pi a}{h} < u < \frac\pi2$: $A_{a,h-a}(2iu) < 0$ (the Type-B breather regime of Harder, Iskandar and McGhee `hep-th/9409035`, eq. (3.16)). Here $\xi_2 = \bar\xi_1$ has $A_{a,h-a}(2iu)E_1 E_2 < 0$, so $\tau_j$ is real and changes sign (singular along curves), whereas $\xi_2 = \bar\xi_1 + i\pi$ has $A_{a,h-a}(2iu)E_1 E_2 > 0$ and is regular.
       - At $u = \frac\pi2 - \frac{\pi a}{h}$: $A_{a,h-a}(2iu) = 0$, the $E_1 E_2$ term vanishes, and the solution is not a breather of mass $2M_a\cos u$.
       - For $0 < u < \frac\pi2 - \frac{\pi a}{h}$: $A_{a,h-a}(2iu) > 0$! Consequently, for $\xi_2 = \bar\xi_1$ (`exr-breather-real`, `11-reality.qmd:87`), $A_{a,h-a}(2iu)\lvert E_1\rvert^2 > 0$, and for $0 < u < \arccos\big(\sqrt{2}\sin\frac{\pi a}{h}\big)$ the quadratic $\tau_j = 1 + 2\operatorname{Re}(\omega^{ja}E_1) + A_{a,h-a}(2iu)\lvert E_1\rvert^2$ has negative discriminant $4 - 4A_{a,h-a}(2iu) < 0$ and stays **strictly positive** ($\tau_j > 0$ everywhere on $\mathbb R^2$) without changing sign! Thus `@exr-breather-real` (`11-reality.qmd:87`) and the caption of `@fig-breather-theta` (`11-reality.qmd:28`) are **wrong** when $a \ne h/2$ unless $u$ is restricted to $(\frac\pi2 - \frac{\pi a}{h}, \frac\pi2)$ where $A_{a,h-a}(2iu) < 0$ (or $a = h/2$).

3. **Soliton–breather scattering (`11-reality.qmd:62`), excited solitons (`10-solitons.qmd:138`), and the missing-charge problem (`11-reality.qmd:64–75`):**
   - Verified in Phase 1 Rank 20 against Takács and Watts (`hep-th/9810006`, Sections 3.5, 4.3–4.4, and 5), `part5-semiclassics/17-thimbles.qmd` (`#sec-which-saddles-contribute`, `@cnj-missing-charge`), and `part9-foundations/44-twelve-gaps.qmd` (`#sec-critique`, `#sec-sector-status`): all citations and cross-references are accurate.

### Verdict

- **Verdict:** **wrong** on two counts (found in Phase 1 Ranks 3 and 7):
  1. `11-reality.qmd:60–61` and `11-reality.qmd:90–92` (`#exr-singular-points`): need the qualifications that real-rapidity 1-soliton zeros are lines (not isolated points), that local charges are unaffected only via contour deformation around the double pole, and that the topological charge jumps by $\pm\vec\alpha_j^\vee$ when a zero crosses the real axis.
  2. `11-reality.qmd:28` (`#fig-breather-theta`) and `11-reality.qmd:87` (`#exr-breather-real`): for a neutral $(a, h-a)$ breather in $a_n^{(1)}$, $A_{a,h-a}(2iu) < 0$ holds on $(\frac\pi2 - \frac{\pi a}{h}, \frac\pi2)$ (which is all of $(0,\frac\pi2)$ only when $a=h/2$); in `@exr-breather-real`, $\tau_j$ changes sign for $\xi_2 = \bar\xi_1$ only when $A_{a,h-a}(2iu) < 0$.

### Paste-ready replacement text

- For `11-reality.qmd:60–61` and `11-reality.qmd:90–92` (`#exr-singular-points`), see Phase 1 Rank 7 above:
  - In `part3-classical/11-reality.qmd`, lines 60–61, replace:
    ```markdown
    Zeros of $\tau_j$ in the real $(x,t)$ plane are generically isolated points: one complex condition, two real conditions. Non-generic data, such as a single soliton with tuned $\operatorname{Im}\xi$ or a breather with $\xi_2=\bar\xi_1$, give curves of zeros. There the energy density has double poles with zero residue, so integrated charges are unaffected, and the Hirota data continue analytically through the singularity.
    ```
    with:
    ```markdown
    Zeros of $\tau_j$ in the real $(x,t)$ plane are generically isolated points (one complex condition, two real conditions, when $x$ and $t$ enter with linearly independent complex coefficients, as at complex rapidity or in multi-soliton collisions). Non-generic data, such as a real-rapidity single soliton with tuned $\operatorname{Im}\xi$ or a breather with $\xi_2=\bar\xi_1$ and $A_{a\bar a}(2iu)<0$, give curves of zeros. At a zero of $\tau_j$, the energy density has a double pole in $x$ with zero residue, so local conserved charges defined by deforming the integration contour around the pole are unaffected, and the Hirota data continue analytically through the singularity (whereas the topological charge defined along the real axis jumps by a coroot $\pm\vec\alpha_j^\vee$ when a zero crosses the real line).
    ```
  - In `part3-classical/11-reality.qmd`, lines 90–92 (`#exr-singular-points`), replace:
    ```markdown
    ::: {#exr-singular-points}
    Show that a zero of $\tau_j$ in the real $(x,t)$ plane is generically an isolated point, and that near it the energy density has a double pole with zero residue, so that integrated charges are unaffected.
    :::
    ```
    with:
    ```markdown
    ::: {#exr-singular-points}
    Show that when $\partial_x\tau_j$ and $\partial_t\tau_j$ are linearly independent over $\mathbb R$, a zero of $\tau_j$ in the real $(x,t)$ plane is an isolated point, and that near any zero of $\tau_j$ the energy density $T_{00}=-\partial_x^2C$ has a double pole in $x$ with zero residue, so that the energy integrated along a contour deformed around the pole is unaffected.
    :::
    ```
- In `part3-classical/11-reality.qmd`, line 28 (`#fig-breather-theta`) and lines 86–88 (`#exr-breather-real`), replace:
  ```markdown
  ![Left: $\Theta$ reflects the rapidity in the real axis and exchanges $a$ and $\bar a$, so the constituents $(a,\theta_0+iu)$ and $(\bar a,\theta_0-iu)$ of a breather are exchanged. Right: the classical breather mass $2M_a\cos u$, for $u\in(0,\pi/2)$, fills a continuous band. `textbook_code/figures/ch11_f16_breather_theta.py` checks that the breather energy and momentum are real, that $E^2-P^2=(2M_a\cos u)^2$ and that $\Theta$ maps the pair to itself.](figures/ch11_f16_breather_theta.svg){#fig-breather-theta fig-alt="Left: two points in the complex rapidity plane at the same real part and opposite imaginary parts, joined by a double arrow. Right: a cosine curve falling from 2 to 0."}
  ```
  with:
  ```markdown
  ![Left: $\Theta$ reflects the rapidity in the real axis and exchanges $a$ and $\bar a$, so the constituents $(a,\theta_0+iu)$ and $(\bar a,\theta_0-iu)$ of a breather are exchanged. Right: the classical breather mass $2M_a\cos u$, for $u\in(0,\pi/2)$, fills a continuous band (for a regular neutral $(a,\bar a)$ breather in $a_n^{(1)}$ with $a\le h/2$, $A_{a\bar a}(2iu)<0$ restricts $u$ to $(\frac\pi2-\frac{\pi a}h,\frac\pi2)$ [@harder1995]). `textbook_code/figures/ch11_f16_breather_theta.py` checks that the breather energy and momentum are real, that $E^2-P^2=(2M_a\cos u)^2$ and that $\Theta$ maps the pair to itself.](figures/ch11_f16_breather_theta.svg){#fig-breather-theta fig-alt="Left: two points in the complex rapidity plane at the same real part and opposite imaginary parts, joined by a double arrow. Right: a cosine curve falling from 2 to 0."}
  ```
  and in `@exr-breather-real` (`11-reality.qmd:86–88`), replace:
  ```markdown
  ::: {#exr-breather-real}
  Show that the asymptotic data of a breather are $\Theta$-invariant. If moreover $\xi_2=\bar\xi_1$, show that $\bar\tau_j=\tau_j$, so that the field is purely imaginary where all $\tau_j>0$, and that the $\tau_j$ change sign, so this solution is singular along curves in the $(x,t)$ plane. For sine-Gordon, show that the real breather (-@eq-sg-breather) instead has $\xi_2=\bar\xi_1+i\pi$. (A numerical check is in `textbook_code/toda_classical.py`.)
  :::
  ```
  with:
  ```markdown
  ::: {#exr-breather-real}
  Show that the asymptotic data of a breather are $\Theta$-invariant. If moreover $\xi_2=\bar\xi_1$, show that $\bar\tau_j=\tau_j$, so that the field is purely imaginary where all $\tau_j>0$, and that in the breather regime $A_{a\bar a}(2iu)<0$ ($\frac\pi2-\frac{\pi a}h<u<\frac\pi2$ for $a\le h/2$) the $\tau_j$ change sign, so this solution is singular along curves in the $(x,t)$ plane. For sine-Gordon, show that the real breather (-@eq-sg-breather) instead has $\xi_2=\bar\xi_1+i\pi$. (A numerical check is in `textbook_code/toda_classical.py`.)
  :::
  ```

---

## Item 9. Exercises in chapters 9–11

### Overlap

- **Overlap with Phase 1:** **in both**.
  - All 17 exercises across Chapters 9–11 (`6` in Chapter 9, `7` in Chapter 10, `4` in Chapter 11) were checked in Phase 1 (across Ranks 3, 5, 7, 8, 9, 11, 12, 14, 16, 19 and the scripts `check_ch09_lagrangian.py`, `check_ch10_solitons.py`, `check_ch11_reality_breathers.py`).

### Check

Summary of independent solutions and verdicts for all 17 exercises:

1. **`#exr-lax` (`09-lagrangian.qmd:164–166`):** **correct.** Solved in `check_ch09_lagrangian.py` Check 2: with $\vec\alpha_1 = -\vec\alpha_0 = \sqrt{2}$, $c_0 = c_1 = 1$, $2\times 2$ matrices give $[\partial_+ + \mathcal A_+, \partial_- + \mathcal A_-] = -\frac{1}{\sqrt{2}}\big(\beta\partial_+\partial_-\phi + 2\sqrt{2}m^2\sinh(\sqrt{2}\beta\phi)\big)\sigma_3 = 0$.
2. **`#exr-spin3` (`09-lagrangian.qmd:168–170`):** **correct** as an exercise prompt, once line 96 (`09-lagrangian.qmd:96`) is fixed so the student replaces $\beta\to i\beta$ ($\beta^2\to -\beta^2$ and $\cosh\beta\phi\to\cos\beta\phi$). On the sine-Gordon kink (`eq-sg-kink`), $Q_{\mp 3} = -\frac{16m^3}{3\beta^2}e^{\mp 3\theta}$ (or $-\frac{32m^3}{3\beta^2}e^{\mp 3\theta}$ with $\mathcal X_2 = m^2(\cos\beta\phi - 1)(\partial_\pm\phi)^2$).
3. **`#exr-dn-masses` (`09-lagrangian.qmd:172–174`):** **correct.** Using $\sum_{a=1}^{n-2}\sin^2\frac{\pi a}{2n-2} = \frac{n-2}{2}$, we get $\sum_{a=1}^{n-2} 8m^2\sin^2\frac{\pi a}{2n-2} + 2m^2 + 2m^2 = (4n-4)m^2 = 2hm^2$, and for $d_4^{(1)}$ ($h=6$) we get $m_1 = m_s = m_{s'} = \sqrt{2}m$ and $m_2 = \sqrt{6}m$, matching `@exr-d4-masses` (`05-root-systems.qmd:113`).
4. **`#exr-cn-masses` (`09-lagrangian.qmd:176–178`):** **correct.** Folding $a_{2n-1}^{(1)}$ ($h=2n$) under $a\leftrightarrow 2n-a$ retains the symmetric combinations of $(a, 2n-a)$ ($a=1,\dots,n-1$) and species $n$, with masses $m_a = 2m\sin\frac{\pi a}{2n}$ ($a=1,\dots,n$).
5. **`#exr-shg-cubic` (`09-lagrangian.qmd:180–182`):** **correct.** Solved in `check_ch09_lagrangian.py` Check 4: for $a_1^{(1)}$, $n_0\alpha_0^3 + n_1\alpha_1^3 = (\sqrt{2})^3 + (-\sqrt{2})^3 = 0$; for $a_2^{(1)}$ in the $\mathbb Z_3$ charge basis, $C_{112} = C_{122} = 0$ and $\lvert C_{111}\rvert = \lvert C_{222}\rvert = 3m^2\beta = \frac{4\beta}{\sqrt{3}}\Delta_{111}$ (where $\Delta_{111} = \frac{\sqrt{3}}{4}(\sqrt{3}m)^2 = \frac{3\sqrt{3}}{4}m^2$).
6. **`#exr-dorey-a2` (`09-lagrangian.qmd:184–186`):** **correct.** Solved in `check_ch09_lagrangian.py` Check 4: with node 1 white and node 2 black, $w = s_2 s_1$ has order $3$, $\vec\gamma_1 = \vec\alpha_1$, $\vec\gamma_2 = s_1\vec\alpha_2 = \vec\alpha_1 + \vec\alpha_2$, $\Omega_1 = \{\vec\alpha_1, \vec\alpha_2, -\vec\alpha_1-\vec\alpha_2\}$, $\Omega_2 = -\Omega_1$, so three roots in $\Omega_1$ (and three in $\Omega_2$) sum to zero while no two roots in $\Omega_1$ sum to a root in $-\Omega_2 = \Omega_1$ with opposite sign unless all three are in the same orbit.
7. **`#exr-bilinear` (`10-solitons.qmd:142–144`):** **correct** (checked in Item 7.1 above).
8. **`#exr-moving-soliton` (`10-solitons.qmd:146–148`):** **correct** (checked in Item 7.1 above).
9. **`#exr-soliton-mass` (`10-solitons.qmd:150–152`):** **correct** (including every step of the hint; checked in Item 7.5 above).
10. **`#exr-charges-a3` (`10-solitons.qmd:154–156`):** **correct** (solved in `check_ch10_solitons.py` Check 1).
11. **`#exr-lie-duality` (`10-solitons.qmd:158–160`):** **correct.** Folding $d_{n+1}^{(1)}$ ($h=2n$, masses $2\sqrt{2}m\sin\frac{\pi a}{2n}$ for $a=1,\dots,n-1$ and $\sqrt{2}m,\sqrt{2}m$ for the two spinors) by the $\mathbb Z_2$ swap of the two spinor nodes keeps the $n-1$ vector masses $2\sqrt{2}m\sin\frac{\pi a}{2n}$ and the symmetric spinor combination of mass $\sqrt{2}m$, which are proportional to $2\sin\frac{\pi a}{2n}$ ($a=1,\dots,n-1$) and $1$, matching the soliton mass ratios of $c_n^{(1)}$.
12. **`#exr-fusion` (`10-solitons.qmd:162–164`):** **wrong** as stated if the student picks $a=b=2$ in $a_2^{(1)}$ ($a+b = 4 > h = 3$), where $m_2 e^{2\pi i/3} + m_2 e^{-2\pi i/3} = -m_1 < 0$; needs the restriction $a+b < h$ (i.e. $a=b=1$ in $a_2^{(1)}$) and the joint limit $\theta_1-\theta_2\to i\pi(a+b)/h, \xi_1,\xi_2\to-\infty$ (see Item 7 above).
13. **`#exr-breather-sg` (`10-solitons.qmd:166–168`):** **correct** (solved in Item 8.2 above).
14. **`#exr-theta-solutions` (`11-reality.qmd:78–80`):** **correct.** For real coupling, $\partial_+\partial_-\vec\phi = -\frac{m^2}{\beta}\sum_j n_j\vec\alpha_j e^{\beta\vec\alpha_j\cdot\vec\phi}$ has real coefficients, so $\vec\phi\mapsto\bar{\vec\phi}$ (optionally composed with a diagram symmetry) maps solutions to solutions, whereas $\vec\phi\mapsto -\bar{\vec\phi}$ flips the sign of $\vec\alpha_j\cdot\vec\phi$ in the exponential and fails unless $-\vec\alpha_j$ is a permutation of the affine simple roots with the same Kac labels (which happens only for $a_1^{(1)}$).
15. **`#exr-theta-interaction` (`11-reality.qmd:82–84`):** **correct** (solved in Item 1.1 above).
16. **`#exr-breather-real` (`11-reality.qmd:86–88`):** **wrong** for $a \ne h/2$ when $0 < u < \frac\pi2 - \frac{\pi a}{h}$, where $A_{a\bar a}(2iu) > 0$ and $\tau_j$ need not change sign; needs the condition $A_{a\bar a}(2iu) < 0$ ($\frac\pi2 - \frac{\pi a}{h} < u < \frac\pi2$ for $a\le h/2$; see Item 8 above).
17. **`#exr-singular-points` (`11-reality.qmd:90–92`):** **wrong** without linear independence of $\partial_x\tau_j, \partial_t\tau_j$ over $\mathbb R$ and contour deformation around the double pole (see Item 8 above).

---

## Overlap table (Phase-2 items 1–9 vs Phase-1 ranks 1–20)

| Phase-2 Item | Phase-1 Rank(s) | Classification | Verdict / Summary of comparison |
|---|---|---|---|
| **Item 1:** $\Theta$ selection rule & Appendix E (`11-reality.qmd`, `E-original-results.qmd:8`) | **Rank 1, Rank 2** | **in both** (prior art & spin-4 check in **Phase 2 only**) | Both found the missing `::: {#prp-...}` statement div/label in Ch. 11 vs Appendix E, the hidden additivity/$\xi_i$-independence and $\operatorname{Re}\cosh\theta_i>0$ assumptions, and the Khastgir–Sasaki negative-energy Cauchy-data distinction. **Phase 2 added** the prior-art check against Fring et al. (`arXiv:1602.05465`, `1608.01691`, `2007.15425`, `2110.06825`) and the spin-4 / general-$s$ $\epsilon_{Q_s}=(-1)^{s-1}$ check. |
| **Item 2:** `@prp-hirota-energy` (`11-reality.qmd:32–41`) | **Rank 1, Rank 4** | **in both** | Both checked `olive1993b`, `zhu1995`, `tw1999`, verified $T_{\mu\nu}^{\rm can}=(\eta_{\mu\nu}\partial^2-\partial_\mu\partial_\nu)C$, and found the missing hypotheses $\operatorname{Re}\cosh\theta_i>0$ and $\prod_{k<l}A_{a_k a_l}\notin\{0,\infty\}$. |
| **Item 3:** Additivity assumption tested (`11-reality.qmd:46`) | **Rank 2** | **in both** | Both checked Freeman (`hep-th/9408092`), constructed the $a_2^{(1)}$ spin-2 current symbolically, and integrated it on 1-solitons and breathers to verify $2iq_1 e^{2\theta_0}\sin 2u$. |
| **Item 4:** Time delays and Wigner sign (`10-solitons.qmd:122`) | **Rank 13** (plus Obs. 1) | **in both** | Both verified $\Delta x, \Delta t$, Wigner's relation, and the sine-Gordon limit; Phase 2 also checked Hollowood's $a_2^{(1)}, a_3^{(1)}$ quantum soliton S-matrix limit and soliton 2's shift. |
| **Item 5:** Lax pair, zero curvature, boost & spins (`09-lagrangian.qmd:32–99`) | **Rank 8, Rank 10, Rank 16** (plus Obs. 2) | **in both** | Both found the two errors at `09-lagrangian.qmd:96` ($\beta^2\to -\beta^2$ without $\cosh\to\cos$, and spin $-3$ vs $+3$ for $\partial_+$) and verified the Lax pair, boost, and exponents mod $kh$. |
| **Item 6:** Masses, couplings, charge basis & Dorey's rule (`09-lagrangian.qmd:102–161`) | **Rank 11, Rank 12, Rank 14** | **in both** (even-$d_{2m}^{(1)}$ error in **Phase 1 only**) | Both checked the mass table, sum rule, area rule, and Dorey's rule; **Phase 1 (Rank 11)** additionally caught that `09-lagrangian.qmd:130` is false for even $d_{2m}^{(1)}$, where the charge eigenbasis of $(s,s')$ *is* real orthonormal. |
| **Item 7:** Solitons, charges, Lie duality & Part VIII pointers (`10-solitons.qmd`) | **Rank 5, Rank 6, Rank 9, Rank 15, Rank 17, Rank 18, Rank 19** | **in both** ($a_{2n}^{(2)}$ degree & fusion errors in **Phase 1 only**) | Both checked `eq-hirota-bilinear`, `eq-an-soliton`, embedded kinks, $c_2^{(1)}$ $\tau_j$, topological charges, `@exr-soliton-mass`, and Lie duality. **Phase 1** additionally caught the $a_{2n}^{(2)}$ exception to $\lvert O\rvert n_j^\vee$ (`Rank 6`), the $\tau_1=\tau_3$ parent-index slip in $c_2^{(1)}$ (`Rank 6`), and the missing $a+b<h$ condition and joint-limit issue in soliton fusion (`Rank 19`), plus checked `eq-q-physical` (`Rank 18`). |
| **Item 8:** Singular solutions, breathers & missing charges (`10-solitons.qmd:128–139`, `11-reality.qmd:58–75`) | **Rank 3, Rank 7, Rank 20** | **in both** (breather $u$-range & topological-charge jump in **Phase 1 only**) | Both checked double poles with zero residue, `tw1999`, and missing charges. **Phase 1** additionally caught that $A_{a,h-a}(2iu)$ changes sign at $u=\frac\pi2-\frac{\pi a}{h}$ when $a\ne h/2$ (`Rank 3`), that real-rapidity 1-soliton zeros are lines rather than isolated points (`Rank 7`), and that the topological charge jumps by $\pm\vec\alpha_j^\vee$ when a zero crosses the real axis (`Rank 7`). |
| **Item 9:** Exercises in chapters 9–11 | **Ranks 3, 5, 7, 8, 9, 11, 12, 14, 16, 19** | **in both** | All 17 exercises were checked in Phase 1 and Phase 2. |

---

## What each list found

### (a) Genuine errors found in Phase 1 that were not targeted in the Phase-2 item prompts

1. **`09-lagrangian.qmd:96` (`#eq-spin3-current`, Phase 1 Rank 8):** The text states *"For sine-Gordon, replace $\beta^2$ by $-\beta^2$"*, which leaves $\cosh(\beta\phi)$ instead of $\cos(\beta\phi)$ and violates conservation on the sine-Gordon equation; moreover, because $\partial_+ = \partial_t+\partial_x$ multiplies a mode of rapidity $\theta$ by $e^{-\theta}$, the charge built from `eq-spin3-current` has spin $-3$ and scales as $e^{-3\theta}$ rather than $e^{+3\theta}$. (Item 5 asked only to verify `eq-spin3-current` for sinh-Gordon and evaluate $Q_3$ on the kink.)
2. **`09-lagrangian.qmd:130` (`#sec-masses-couplings`, Phase 1 Rank 11):** The text claims that for degenerate masses in $(a_n, d_n, e_6)$ the charge eigenvectors are complex and that *"In a real orthonormal basis of a degenerate pair neither (-@eq-area-rule) nor @prp-dorey-rule holds"*. This is false for $d_{2m}^{(1)}$ with even $n=2m$ (including $d_4^{(1)}$, plotted in `@fig-dorey-d4` on the same page), where the $\mathbb Z_2\times\mathbb Z_2$ centre characters are real ($\pm 1$), $\bar a = a$, and the charge eigenbasis of $(s,s')$ *is* a real orthonormal basis in which both `eq-area-rule` and `@prp-dorey-rule` hold. (Item 6 asked only to verify that a real orthonormal basis fails for $a_n^{(1)}$ and $e_6^{(1)}$.)
3. **`10-solitons.qmd:128–136`, `11-reality.qmd:28` (`#fig-breather-theta`), and `11-reality.qmd:87` (`#exr-breather-real`) (Phase 1 Rank 3):** For a neutral $(a, h-a)$ breather of $a_n^{(1)}$ with $a < h/2$, $A_{a,h-a}(2iu) = 1 - \sin^2(\pi a/h)/\cos^2 u$ vanishes at $u = \frac\pi2 - \frac{\pi a}{h} \in (0,\frac\pi2)$ and is positive for $0 < u < \frac\pi2 - \frac{\pi a}{h}$. Thus a regular neutral breather (and the sign change of $\tau_j$ at $\xi_2 = \bar\xi_1$ claimed in `@exr-breather-real`) holds only on the sub-interval $\frac\pi2 - \frac{\pi a}{h} < u < \frac\pi2$, not on all of $(0,\frac\pi2)$.
4. **`10-solitons.qmd:76–77` (`#sec-vertex-solitons`, Phase 1 Rank 6):** The polynomial degree formula $\lvert O\rvert n_j^\vee$ fails for the self-dual twisted algebra $a_{2n}^{(2)}$ folded from $a_{2n}^{(1)}$, where the degree is $\lvert O\rvert n_j^{\rm parent} = 2$ at every node while $\lvert O\rvert n_j^\vee$ gives $4$ at nodes $0,\dots,n-1$; and the $c_2^{(1)}$ example writes $\tau_1 = \tau_3$ using the parent $a_3^{(1)}$ node index $3$ without explanation.
5. **`10-solitons.qmd:126` and `10-solitons.qmd:162–164` (`#exr-fusion`, Phase 1 Rank 19):** The fusion identity $m_a e^{i\pi b/h} + m_b e^{-i\pi a/h} = m_{a+b}$ requires $a+b < h$ (for $a+b > h$, $\sin\frac{\pi(a+b)}{h} < 0$ and the fused species $a+b-h$ arises at $\theta_1 = \theta - i\pi(h-b)/h, \theta_2 = \theta + i\pi(h-a)/h$), and requires a joint limit $\theta_1-\theta_2 \to i\pi(a+b)/h$ with $\xi_1,\xi_2\to -\infty$.
6. **`11-reality.qmd:60–61` and `11-reality.qmd:90–92` (`#exr-singular-points`, Phase 1 Rank 7):** A single soliton with real rapidity $\theta$ depends only on $x\cosh\theta - t\sinh\theta$, so its real zeros are lines, never isolated points in $\mathbb R^2$; local conserved charges are unaffected by a zero of $\tau_j$ on $\mathbb R$ only when the integration contour is deformed into $\mathbb C$ around the double pole; and the topological charge jumps by $\pm\vec\alpha_j^\vee$ when a zero crosses the real axis.
7. **`11-reality.qmd:43–45` vs `appendices/E-original-results.qmd:8` (Phase 1 Rank 2):** The $\Theta$ selection rule in `11-reality.qmd` is an unlabelled bullet point without a `::: {#prp-...}` statement div or status label, whereas `appendices/E-original-results.qmd:8` lists it with the status label `**[Theorem, given additivity of the higher charges]**` and cites `@sec-reality-energy-momentum` instead of a statement id.
8. **`11-reality.qmd:47` and `10-solitons.qmd:180` (Phase 1 Rank 1):** Khastgir and Sasaki's (`hep-th/9507001`) smooth negative-energy solutions are Cauchy configurations with $\vec\phi(x,0)=0$ and purely imaginary $\partial_t\vec\phi(x,0) = i\vec f(x)$, which are strictly $\Theta$-invariant ($\Theta\vec\phi = \vec\phi$); they are excluded by restricting to the Hirota multi-soliton sector with $\operatorname{Re}\cosh\theta_i > 0$, not by $\Theta$-invariance alone.

### (b) Genuine errors found through Phase-2 items 1–9 that Phase 1 missed

1. **Missing prior-art citations and overbroad originality claim for the $\Theta$ selection rule (`11-reality.qmd:96` and `appendices/E-original-results.qmd:8`, found in Item 1):**
   Phase 1 did not search the external post-2015 literature of Fring and collaborators (`arXiv:1602.05465`, `arXiv:1608.01691`, `arXiv:2007.15425`, `arXiv:2110.06825`). Reading those four papers in Phase 2 revealed that Cen–Fring (`arXiv:1602.05465`) and Cen–Correa–Fring (`arXiv:1608.01691`) already introduced the two-step mechanism combining antilinear $\mathcal{PT}$ symmetry ($\phi\mapsto -\bar\phi, (x,t)\mapsto(-x,-t)$) on single solitons with asymptotic additivity over constituent solitons to explain reality of energy and higher conserved charges for multi-soliton configurations at real rapidities, and Correa–Fring–Taira (`arXiv:2110.06825`) applied this to Bullough–Dodd ($a_2^{(2)}$) solitons and breathers. Therefore `11-reality.qmd:96` (`"The formulation as a selection rule from $\Theta$ is new here"`) and `appendices/E-original-results.qmd:8` must cite those papers and state specifically what is new in Chapter 11 (namely, formulating the selection rule on multicomponent Hirota data $\{(a_i,\theta_i)\}$ under $\Theta:(a,\theta)\mapsto(\bar a,\bar\theta)$, including complex-rapidity conjugate pairs/breathers and the charge-conjugation sign $\epsilon_Q = (-1)^{s-1}$ for $C$-odd higher charges).

### (c) What each approach was better at finding

- **Phase 1 (unprompted adversarial risk-ranking from the chapter text and local repository):**
  Was substantially better at finding **mathematical boundary-case failures, quantifier/domain errors, convention/sign slips, and internal cross-chapter inconsistencies**. Because Phase 1 systematically questioned every universal quantifier in Chapters 9–11 ("every degenerate pair", "every folded theory", "$u\in(0,\pi/2)$", "every $a,b$", "generically isolated points", "integrated charges are unaffected", "replace $\beta^2$ by $-\beta^2$"), it caught 8 mathematical and structural errors across the three chapters that a guided checklist confirming the main derivations did not flag.
- **Phase 2 (brief-guided checklist with targeted external prior-art pointers):**
  Was essential for **targeted external prior-art verification** — specifically, pointing to the four papers of Fring et al. (`arXiv:1602.05465`, `1608.01691`, `2007.15425`, `2110.06825`) that were not cited in `references.bib` or stored in `/home/gustav/referee-copy/papers/`. It also prompted valuable extensions of the computational checks (such as the spin-4 parity check in $a_2^{(1)}$, the $3\times 3$ Bullough–Dodd Lax pair check, and the $a_4^{(1)}, e_6^{(1)}$ coupling checks).

---

## Claims that could not be checked, and why

None. Every item and sub-question in `briefs/referee_classical_09_11_items.md` was checked in full by independent derivation, symbolic/numerical scripts in `referee_scratch/classical_09_11/`, and primary source reading.

---

## Other observations (unverified)

1. **`part3-classical/10-solitons.qmd:122` (`#sec-multi-solitons`): Overload of the symbol $E_1$ within a single paragraph.**
   In the first sentence of the "Time delays" paragraph, $E_1$ and $E_2$ denote the Hirota exponentials of solitons 1 and 2 (`"either $E_2\to0$, and $\tau_j\simeq1+\omega^{ja_1}E_1$, or $E_2\to\infty$"`). Five sentences later in the same paragraph, $E_1$ denotes the energy of soliton 1 (`"the shift is $E_1\Delta x=-p_1\Delta t=-(2h/\beta^2)\ln A_{a_1a_2}$, with $E_1,p_1$ the energy and momentum of soliton 1"`), despite the notation note at line 61 (`"In this chapter and the next... $E$ is the exponential in (-@eq-an-soliton), not an energy"`). Writing $\mathcal E_1$ (or $M_{a_1}\cosh\theta_1\,\Delta x$) for the energy of soliton 1 at line 122 would avoid the collision.
2. **`part3-classical/09-lagrangian.qmd:90` (`#eq-spin3-current`) vs `09-lagrangian.qmd:165` (`#exr-lax`): Different normalizations of $\beta$ and $m$ for sinh-Gordon within the same chapter.**
   At line 90, sinh-Gordon is written with $V = \frac{m^2}{\beta^2}(\cosh\beta\phi - 1)$ and $\partial_+\partial_-\phi = -\frac{m^2}{\beta}\sinh\beta\phi$ (where the particle mass is $m$), whereas in `#eq-imaginary-potential` (line 24) and `@exr-lax` (line 165), $a_1^{(1)}$ has $\vec\alpha_1 = -\vec\alpha_0 = \sqrt{2}$, $V_{\rm real} = \frac{2m^2}{\beta^2}(\cosh(\sqrt{2}\beta\phi)-1)$, and particle mass $m_1 = 2m$. Adding a parenthetical `(with $\sqrt2\beta$ and $2m$ written as $\beta$ and $m$, as in @sec-sg-dictionary)` at line 90 would prevent confusion when comparing `@exr-lax` and `@exr-spin3`.
