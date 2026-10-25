# Referee Report: `lecture-notes/defects/` (Lectures 1–9 and Home-Page Count)

- **Date:** 2026-10-25
- **Model:** Gemini 4 Argon
- **Commit reviewed:** `9b912ee` (`main` in worktree; includes PRs #30, #34, #50, #55; note that PR #56 adding `lecture-notes/defects/09-folding-fusing.qmd` was merged on GitHub after `9b912ee`, so Item 5 evaluates Lecture 9's mathematical claims from `briefs/defects_lectures_handoff.md` §7.2, PR #56's specification, and `part8-defects/43-defects-and-boundaries.qmd`, while flagging that `09-folding-fusing.qmd` is not yet in commit `9b912ee`)
- **Independent verification scripts:**
  - `referee_scratch/lectures_defects/check_l1_l5.py` (Lectures 1–5)
  - `referee_scratch/lectures_defects/check_l6.py` (Lecture 6)
  - `referee_scratch/lectures_defects/check_l7_l9.py` (Lectures 7–9)

---

## Summary Table of Verdicts

| Item | Lecture / Statement | Current Label | Verdict | `[Theorem]` Label Justified? | Notes / Required Action |
|---|---|---|---|---|---|
| **1.1** | L1: `thm-l1-nogo` (`01-mirrors-and-windows.qmd:241`) | `[Theorem]` | **correct** (minor proof gap on right-incidence) | Yes (once 1 sentence citing defect unitarity is added) | Proof establishes $R_a \equiv 0$ or $\mathbb{T}_a \equiv 0$; add 1 sentence using (-@eq-l1-unitarity) to deduce $\bar{R}_a \equiv 0$ or $\bar{\mathbb{T}}_a \equiv 0$ and exclude $R_a = \mathbb{T}_a = 0$. |
| **1.1b** | L1: DMS (`hep-th/9409076`) eqs. (4.6), (4.8) note (`01-mirrors-and-windows.qmd:80`) | Note | **correct** | N/A | Independently verified that $\mathcal{L}_{\rm def} = -\frac{\lambda}{2}\delta(x)\phi^2$ gives $\tan\chi = \frac{\lambda}{2m} = \frac{g}{2m}$, so DMS (4.6) and (4.8) have a factor-of-2 misprint $g/(4m)$. |
| **1.2** | L2: `prp-l2-bianchi` (`02-backlund-machine.qmd:245`) & geometry | `[Established]` | **correct** | N/A (`[Established]`) | Necessity proved in text; sufficiency cited to Rogers–Schief. Beltrami pseudosphere, two-soliton, and breather formulas verified. |
| **1.3a** | L3: `thm-l3-momentum` (`03-freezing.qmd:241`) | `[Theorem]` | **correct** | Yes | Proof is complete and verified symbolically; add 1 sentence at `03-freezing.qmd:183` noting why a general linear velocity term $A\Phi_t + B\Psi_t$ reduces to $k(\Phi,\Psi)\Psi_t$ modulo a total time derivative. |
| **1.3b** | L3: `prp-l3-which` (`03-freezing.qmd:309`) | `[Theorem]` | **wrong** (minor overclaim in $\kappa = 0$ case) | Yes (with 1-line fix to $\kappa = 0$ case) | When $\kappa = 0$ ($F'' = G'' = 0$), $F(\Sigma) = a_1\Sigma + a_0$ and $G(\Delta) = b_1\Delta + b_0$ are affine-linear: if $a_1 = 0, a_0 b_1 \neq 0$, $V_L(\Phi) = 2a_0 b_1\Phi + c$ is a non-zero linear potential ($\square\Phi = -2a_0 b_1 \neq 0$) while $\mathcal{B} = a_0(\Phi+\Psi) + \frac{1}{2}b_1(\Phi-\Psi)^2$ depends on both $\Phi+\Psi$ and $\Phi-\Psi$. Fix statement and proof of the $\kappa = 0$ case. |
| **1.4** | L4: `prp-l4-delay` (`04-what-gets-through.qmd:127`), linear waves, half-solitons | `[Theorem]` | **correct** | Yes | All formulas ($z(\theta)$, $\mathbb{T}(\theta) = z(\theta+i\pi/2)$, $z_{\rm odd}^2 = A$, $\mathbb{T}_{\rm odd}^2 = f_{1/2}$, energy identities) verified symbolically. |
| **1.5** | L5: `prp-l5-cycle` (`05-why-the-cycle.qmd:221`) & `book:` links in L3, L5 | `[Theorem]` | **correct** | Yes | Proof of $A_{\alpha\beta} = J_{\alpha\beta} - C_{\alpha\beta}$ and type-I condition $C_{\alpha\beta} = 0$ verified symbolically and for $d_4^{(1)}$. Update `book:` links in `03-freezing.qmd` and `05-why-the-cycle.qmd` to `book:prp-dc-types`, `book:eq-dc-halves`, `book:eq-dc-lightcone`, `book:eq-dc-components`. |
| **2.1** | L6: `thm-l6-gauge` (`06-classical-algebra.qmd:123`) | `[Theorem]` | **correct** | Yes | Complete 4-step proof; verified symbolically. |
| **2.2** | L6: `thm-l6-frozen` (`06-classical-algebra.qmd:210`) | `[Theorem]` | **correct** | Yes | Complete proof for both type-I$_-$ and type-I$_+$; verified symbolically. |
| **2.3** | L6: `thm-l6-monodromy` (`06-classical-algebra.qmd:273`) | `[Theorem]` | **correct** | Yes | Complete proof from transport across left/right halves and gauge relation. |
| **2.4** | L6: `prp-l6-bulk-r` (`06-classical-algebra.qmd:404`) & `thm-l6-sklyanin` (`:424`) | `[Established]` / `[Theorem]` | **correct** | Yes | All 16 components of $\{L_d(\theta_1)\stackrel{\otimes}{,}L_d(\theta_2)\} = [r(\theta_1-\theta_2), L_d(\theta_1)\otimes L_d(\theta_2)]$ verified symbolically. |
| **2.5** | L6: Spin-$\mp 3$ defect terms $\mathcal{D}_{\mp 3}$ (`06-classical-algebra.qmd:368`) & $a_r^{(1)}$ $L_d$ (`:518`) | Equations | **correct** | N/A | Exact differential identity $\partial_t\mathcal{D}_{-3} = \mathcal{T}_{-4}(\Phi) - \mathcal{T}_{-4}(\Psi)$ verified symbolically under (-@eq-l6-defect), along with $a_r^{(1)}$ Lax gauge equation. |
| **3.1** | L7: `thm-l7-half` (`07-quantizing.qmd:117`) | `[Theorem]` | **correct** | Yes | Complete algebraic proof that $[e_1, f_0]=0$ forces $c_n$ constant, contradicting $[e_1, f_1] = (k_1-k_1^{-1})/(q-q^{-1})$. |
| **3.2** | L7: `thm-l7-T` (`07-quantizing.qmd:219`) & `prp-l7-triangle` (`:263`) | `[Theorem]` / `[Established]` | **correct** | Yes | All 8 linear intertwining equations solved and verified symbolically; `prp-l7-triangle` properly labelled `[Established]`. |
| **3.3** | L7: `thm-l7-crossing` (`07-quantizing.qmd:386`) & `thm-l7-poles` (`:443`) | `[Theorem]` | **correct** | Yes | Crossing sign condition $q = -e^{-i\pi\lambda}$, integral representation (-@eq-l7-Phi), poles, residues, and Wigner delay slope $\delta'(0) = \frac{\lambda}{\pi}(\ln\frac{8\lambda}{\pi}+\gamma_E)+o(1)$ all verified. |
| **3.4** | L7: `thm-l7-breather` (`07-quantizing.qmd:572`) | `[Theorem]` | **correct** | Yes | Complete proof via crossing and Fourier-kernel cancellation of $2\sinh\frac{\pi k}{2\lambda}$; verified to $4.58\times 10^{-16}$. |
| **4.1** | L8: `thm-l8-turned` (`08-q-operators.qmd:56`) | `[Theorem]` | **correct** | Yes | Complete proof via Euclidean continuation $t = -i\tau$ and swap $(\tau, x) \mapsto (x, t)$. |
| **4.2** | L8: `thm-l8-commute` (`08-q-operators.qmd:181`) & `prp-l8-dd` (`:192`) | `[Theorem]` / `[Numerical]` | **correct** (minor technical note on $X^{-1}$) | Yes (with 1-phrase note on $e^{2\lambda(\theta-\theta_d)} \neq -1$) | Proof of `thm-l8-commute` uses $X^{-1}$ where $X = \mathbb{L}_{a\mathcal{W}}(\theta-\theta_d)$, which is invertible for $e^{2\lambda(\theta-\theta_d)} \neq -1$ and extends to $e^{2\lambda(\theta-\theta_d)} = -1$ by analyticity/continuity in $\theta$. |
| **4.3** | L8: `thm-l8-tq` (`08-q-operators.qmd:239`) & `prp-l8-baxter` (`:284`) | `[Theorem]` | **correct** | Yes | Short exact sequence on $\mathcal{S} = \ker X$ and quotient verified symbolically and numerically for $N=2,3$ across all charge sectors to $2.9\times 10^{-15}$. |
| **5.1** | L9: `09-folding-fusing.qmd` (PR #56) & Ch. 43 (`43-defects-and-boundaries.qmd`) | `[Theorem]` / `[Proposition]` | **correct** (file not yet in `9b912ee`) | Yes (for `thm-l9-dirichlet` integral proof) | Trigonometric identity behind the integral proof of `thm-l9-dirichlet` (`prp-defect-fusion-k` (3)) and second Borel half $\kappa = CB^{-1} = q^{-2n-2-s}$ verified. See Item 5 for `43-defects-and-boundaries.qmd` chirality check. |
| **5.2** | Home page & index (`index.qmd:38,46,49`, `lecture-notes/defects/index.qmd:21`) | Prose / Table | **wrong** (outdated counts/status) | N/A | `index.qmd` still says "Lectures 1–6... written" and "53 so far". Actual count is **73 exercises** in Lectures 1–8 (at `9b912ee`, or **83** once Lecture 9 is merged) and **14 guided problems** in `defects-gemini`. Exact replacement text given below. |

---

## Item 1: Lectures 1–5 (`01-mirrors-and-windows.qmd` to `05-why-the-cycle.qmd`, PR #30)

### 1.1 Lecture 1 (`01-mirrors-and-windows.qmd`)
- **`thm-l1-nogo` (`01-mirrors-and-windows.qmd:241`, `[Theorem]`): Verdict: `correct` (with minor proof completion for right-incidence amplitudes).**
  - **Verification:** In a diagonal bulk theory with Zamolodchikov–Faddeev relations $A_a(\theta_1)A_b(\theta_2) = S_{ab}(\theta_1-\theta_2)A_b(\theta_2)A_a(\theta_1)$ and $A_a(\theta)D = R_a(\theta)A_a(-\theta)D + \mathbb{T}_a(\theta)D A_a(\theta)$, commuting $A_a(\theta_1)A_b(\theta_2)D$ in two ways yields four channels ($LL, LR, RL, RR$). The mixed left-right channels $A_a(-\theta_1)D A_b(\theta_2)$ and $A_b(-\theta_2)D A_a(\theta_1)$ give
    $$
    R_a(\theta_1)\mathbb{T}_b(\theta_2)\big[S_{ab}(-\theta_1-\theta_2) - S_{ab}(\theta_1-\theta_2)\big] = 0.
    $$
    Since $S_{ab}(\theta)$ is not constant for any pair $(a,b)$, $S_{ab}(-\theta_1-\theta_2) - S_{ab}(\theta_1-\theta_2)$ cannot vanish on an open set, forcing $R_a(\theta_1)\mathbb{T}_b(\theta_2) \equiv 0$ for all $a,b$. This proves that either $R_a \equiv 0$ for all $a$ or $\mathbb{T}_a \equiv 0$ for all $a$.
  - **Minor gap to close:** The theorem statement (`01-mirrors-and-windows.qmd:243`) concludes that "Every integrable defect is either purely reflecting or purely transmitting", which also requires knowing that (i) $R_a$ and $\mathbb{T}_a$ do not *both* vanish, and (ii) the right-incident amplitudes $\bar{R}_a, \bar{\mathbb{T}}_a$ obey the same alternative ($\bar{R}_a \equiv 0$ when $R_a \equiv 0$, and $\bar{\mathbb{T}}_a \equiv 0$ when $\mathbb{T}_a \equiv 0$). Both follow immediately from defect unitarity (-@eq-l1-unitarity), $R_a(\theta)R_a(-\theta) + \mathbb{T}_a(\theta)\bar{\mathbb{T}}_a(-\theta) = 1$ and $R_a(\theta)\mathbb{T}_a(-\theta) + \mathbb{T}_a(\theta)\bar{R}_a(-\theta) = 0$ (together with its right-incident partner $\bar{R}_a(\theta)\bar{R}_a(-\theta) + \bar{\mathbb{T}}_a(\theta)\mathbb{T}_a(-\theta) = 1$).
  - **Replacement text (`lecture-notes/defects/01-mirrors-and-windows.qmd`, lines 256–258):**
    ```markdown
    If some $R_a$ is not identically zero, pick $\theta_1$ with $R_a(\theta_1)\ne0$. For each $b$, the bracket vanishes at most for isolated $\theta_2$, because $S_{ab}$ is not a constant. Hence $\mathbb T_b(\theta_2)=0$ for all $\theta_2$ and all $b$: the defect is a mirror. By the same argument, if some $\mathbb T_b$ is not identically zero, every $R_a$ vanishes: the defect is a window. Defect unitarity (-@eq-l1-unitarity) excludes $R_a=\mathbb T_a=0$ and forces $\bar{\mathbb T}_a=0$ in the first case and $\bar R_a=0$ in the second.
    :::
    ```
- **DMS (`hep-th/9409076`) factor-of-2 observation (`01-mirrors-and-windows.qmd:80`): Verdict: `correct`.**
  - **Verification (`PASS L1.1`, `PASS L1.2` in `check_l1_l5.py`):** Varying $S = \int dt\left[\int_{x\neq 0} dx\,\frac{1}{2}((\partial_\mu\phi)^2 - m^2\phi^2) - \frac{\lambda}{2}\phi(t,0)^2\right]$ gives the jump condition $\phi_x(t,0^+) - \phi_x(t,0^-) = \lambda\phi(t,0)$. For a plane wave of momentum $k = m\sinh\theta$, $R(\theta) = \frac{-i\lambda}{2m\sinh\theta + i\lambda} = \frac{-i\tan\chi}{\sinh\theta + i\tan\chi}$ with $\tan\chi = \frac{\lambda}{2m}$. In Delfino–Mussardo–Simonetti (`hep-th/9409076`), eq. (4.1) writes $\mathcal{L}_{\rm defect} = -\frac{1}{2}g\delta(x)\phi^2$ (so $g = \lambda$), which gives $\tan\chi = \frac{g}{2m}$, whereas their eqs. (4.6) and (4.8) print $\frac{g}{4m}$.

### 1.2 Lecture 2 (`02-backlund-machine.qmd`)
- **`prp-l2-bianchi` (`02-backlund-machine.qmd:245`, `[Established]`): Verdict: `correct`.**
  - **Verification (`PASS L2.1`–`PASS L2.4` in `check_l1_l5.py`):**
    1. Necessity of the algebraic Bianchi relation $\tan\frac{\Phi_3 - \Phi_0}{4} = \frac{a_1+a_2}{a_1-a_2}\tan\frac{\Phi_1-\Phi_2}{4}$ is proved completely in lines 206–234 by subtracting the $\partial_+,\partial_-$ Bäcklund equations around the plaquette and dividing.
    2. We also verified symbolically that differentiating the algebraic relation $\Phi_3(\Phi_0,\Phi_1,\Phi_2)$ with respect to $\partial_+$ and $\partial_-$ and substituting the first-order Bäcklund relations for $\Phi_1,\Phi_2$ implies that $\Phi_3$ automatically satisfies both Bäcklund equations with $\Phi_1$ (parameter $a_2$) and $\Phi_2$ (parameter $a_1$). Because the lecture text cites Rogers–Schief [@rogers_schief2002] for sufficiency rather than writing out that differentiation, the label `[Established]` is appropriate.
    3. All geometric and multi-soliton formulas in Lecture 2 (`eq-l2-pseudosphere`, `eq-l2-tractrix`, `eq-l2-2soliton`, `eq-l2-breather`) passed symbolic and numerical PDE checks (`max_err = 9.05e-12`).

### 1.3 Lecture 3 (`03-freezing.qmd`)
- **`thm-l3-momentum` (`03-freezing.qmd:241`, `[Theorem]`): Verdict: `correct` (`[Theorem]` justified).**
  - **Verification (`PASS L3.1`, `PASS L3.2` in `check_l1_l5.py`):** In light-cone coordinates $\Sigma = \Phi+\Psi$, $\Delta = \Phi-\Psi$, the momentum flux across $x=0$ for the first-order defect action with $k(\Phi,\Psi) = \Phi$ reduces to $\mathcal{B}_\Sigma\mathcal{B}_\Delta = \frac{1}{2}(V_L(\Phi) - V_R(\Psi))$ at $x=0$, which forces $\mathcal{B}_{\Sigma\Sigma} = \mathcal{B}_{\Delta\Delta} = 0$, hence $\mathcal{B}(\Phi,\Psi) = F(\Phi+\Psi) + G(\Phi-\Psi)$ with $2F'(\Sigma)G'(\Delta) = V_L(\Phi) - V_R(\Psi)$.
  - **Minor exposition improvement (`03-freezing.qmd:183`):** In Step 1 (`03-freezing.qmd:183`), the text writes $L_d = k(\Phi,\Psi)\Psi_t - \mathcal{B}(\Phi,\Psi)$ after noting that $\Phi\Phi_t$ is a total time derivative. More generally, any term linear in both velocities, $A(\Phi,\Psi)\Phi_t + B(\Phi,\Psi)\Psi_t$, differs from $k(\Phi,\Psi)\Psi_t$ by a total time derivative $\frac{d}{dt}S(\Phi,\Psi)$ (taking $\partial_\Phi S = A$ and $k = B - \partial_\Psi S$). Adding one sentence at line 183 makes the generality of Step 1 explicit:
    ```markdown
    **Step 1: no time derivatives of the jump.** $P_{\text{bulk}}$ contains $\Phi_x\Phi_t$ and $\Psi_x\Psi_t$. For $\dot P_{\text{bulk}}$ to be a total time derivative, the defect conditions at $x=0$ must trade $\Phi_x$ and $\Psi_x$ for $\Phi_t$ and $\Psi_t$. So $L_d$ must be of first order in time, without second derivatives: $L_d=k(\Phi,\Psi)\Psi_t-\mathcal B(\Phi,\Psi)$ (any term $A(\Phi,\Psi)\Phi_t$ can be absorbed into $k(\Phi,\Psi)\Psi_t$ up to a total time derivative $\frac{d}{dt}S(\Phi,\Psi)$ with $S_\Phi=A$, which does not affect the defect conditions).
    ```
- **`prp-l3-which` (`03-freezing.qmd:309`, `[Theorem]`): Verdict: `wrong` (minor overclaim in the $\kappa = 0$ case of the theorem statement and proof).**
  - **What is wrong:** In the proof (`03-freezing.qmd:327`), differentiating $2F'(\Sigma)G'(\Delta) = V_L(\Phi) - V_R(\Psi)$ gives $F'''(\Sigma)/F'(\Sigma) = -G'''(\Delta)/G'(\Delta) = \kappa$. When $\kappa = 0$, $F''$ and $G''$ are constants, so $F(\Sigma)$ and $G(\Delta)$ are quadratic polynomials:
    $$
    F'(\Sigma) = a_1\Sigma + a_0,\qquad G'(\Delta) = b_1\Delta + b_0.
    $$
    Then
    $$
    2F'(\Sigma)G'(\Delta) = \big[2a_1 b_1\Phi^2 + 2(a_1 b_0 + a_0 b_1)\Phi\big] - \big[2a_1 b_1\Psi^2 + 2(a_1 b_0 - a_0 b_1)\Psi\big] + 2a_0 b_0.
    $$
    The theorem statement (`03-freezing.qmd:312–313`) claims that for $\kappa = 0$ the potentials are either "(3) the free massive or tachyonic potential $\frac{1}{2}m^2\Phi^2$, $m^2\neq 0$, with $F$ and $G$ quadratic; and (4) the massless free case $V_L = V_R = 0$, in which $\mathcal{B}$ depends on $\Phi+\Psi$ alone or on $\Phi-\Psi$ alone."
    This misses the case $a_1 = 0, a_0 b_1 \neq 0$ (or $b_1 = 0, a_1 b_0 \neq 0$): for instance, if $F(\Sigma) = a_0(\Phi+\Psi)$ and $G(\Delta) = \frac{1}{2}b_1(\Phi-\Psi)^2$ with $a_0 b_1 \neq 0$, then $V_L(\Phi) = 2a_0 b_1\Phi + c$ and $V_R(\Psi) = 2a_0 b_1\Psi + c$ are *non-zero linear potentials* (massless fields with a constant bulk source, $\Phi_{tt}-\Phi_{xx} = -2a_0 b_1$), while $\mathcal{B}(\Phi,\Psi) = a_0(\Phi+\Psi) + \frac{1}{2}b_1(\Phi-\Psi)^2$ depends on *both* $\Phi+\Psi$ and $\Phi-\Psi$. Even when $a_1 b_1 \neq 0$, $V_L$ and $V_R$ can include linear terms in $\Phi$ and $\Psi$ (which can be absorbed by constant shifts of $\Phi,\Psi$).
  - **Replacement text (`lecture-notes/defects/03-freezing.qmd`, lines 309–329):**
    ```markdown
    ::: {#prp-l3-which}
    ## Which scalar potentials admit a momentum-conserving defect
    **[Theorem]** Up to constant shifts and rescalings of $\Phi$ and $\Psi$ and an overall scale of the potential, the only non-constant potentials $V_L(\Phi)$, $V_R(\Psi)$ that satisfy (-@eq-l3-factor) are: (1) the sine-Gordon and sinh-Gordon potentials, with $V_L$ and $V_R$ equal up to a shift; (2) the Liouville potential $e^{b\Phi}$, joined to another Liouville potential or to the massless free theory $V_R=0$; (3) the quadratic (free massive or tachyonic) potential $\frac12m^2\Phi^2$, $m^2\ne0$, with $F$ and $G$ quadratic; and (4) the linear or massless free case $V''_L=V''_R=0$, where either $V_L(\Phi)=V_R(\Psi)=\mu\Phi$ is linear with one of $F,G$ linear and the other quadratic, or $V_L=V_R=0$ is massless free with $\mathcal B$ depending on $\Phi+\Psi$ alone or on $\Phi-\Psi$ alone.
    :::

    ::: {.proof}
    Apply $\partial_\Sigma\partial_\Delta=\frac14(\partial_\Phi^2-\partial_\Psi^2)$ to (-@eq-l3-factor): $2F''(\Sigma)G''(\Delta)=\frac14\big(V_L''(\Phi)-V_R''(\Psi)\big)$. Apply it once more:

    $$
    2F'''(\Sigma)\,G'''(\Delta)=\frac1{16}\big(V_L^{(4)}(\Phi)-V_R^{(4)}(\Psi)\big).
    $$

    Divide the second equation by the first where $F''G''\ne0$. On the left we have $\big(F'''(\Sigma)/F''(\Sigma)\big)\big(G'''(\Delta)/G''(\Delta)\big)$, a product of a function of $\Sigma$ and a function of $\Delta$. Wait: we do not even need the fourth derivative. Apply $\partial_\Sigma^2-\partial_\Delta^2=\partial_\Phi\partial_\Psi$ directly to (-@eq-l3-factor). The right side of (-@eq-l3-factor) is $V_L(\Phi)-V_R(\Psi)$, a sum of a function of $\Phi$ and a function of $\Psi$, so $\partial_\Phi\partial_\Psi$ kills it:

    $$
    F'''(\Sigma)\,G'(\Delta)-F'(\Sigma)\,G'''(\Delta)=0 .
    $$

    Where $F'$ and $G'$ are non-zero we can divide by $F'(\Sigma)G'(\Delta)$:

    $$
    \frac{F'''(\Sigma)}{F'(\Sigma)}=\frac{G'''(\Delta)}{G'(\Delta)}=\kappa,
    $$

    a separation constant. For $\kappa<0$, $F'$ and $G'$ are combinations of $\sin$ and $\cos$; their product $2F'(\Sigma)G'(\Delta)$, written in $\Phi$ and $\Psi$ by the product-to-sum identities, is a combination of $\cos$ and $\sin$ of $\sqrt{-\kappa}\,\Phi$ and $\sqrt{-\kappa}\,\Psi$, which gives sine-Gordon. For $\kappa>0$, $\sinh$ and $\cosh$ give sinh-Gordon when both exponentials are present in $F'$ and $G'$, and Liouville (on one or both sides) when one exponential is absent. For $\kappa=0$, $F'$ and $G'$ are affine-linear, $F'(\Sigma)=a_1\Sigma+a_0$ and $G'(\Delta)=b_1\Delta+b_0$: if $a_1b_1\ne0$, constant shifts of $\Phi$ and $\Psi$ remove $a_0,b_0$ and give $V_L(\Phi)=V_R(\Psi)=2a_1b_1\Phi^2+c$; if one of $a_1,b_1$ vanishes and the other is non-zero, $V_L$ and $V_R$ are linear with equal slope $2a_0b_1$ or $2a_1b_0$; if $a_1=b_1=0$, one of $a_0,b_0$ must vanish to keep $V_L-V_R$ constant, so $F'=0$ or $G'=0$, which is the massless free case with $\mathcal B$ depending on $\Delta$ alone or on $\Sigma$ alone.
    :::
    ```

### 1.4 Lecture 4 (`04-what-gets-through.qmd`)
- **`prp-l4-delay` (`04-what-gets-through.qmd:127`, `[Theorem]`): Verdict: `correct` (`[Theorem]` justified).**
  - **Verification (`PASS L4.1`, `PASS L4.2` in `check_l1_l5.py`):** Substituting $\tan(\Phi/4) = Y_L$ and $\tan(\Psi/4) = z(\theta)Y_L$ into the frozen Bäcklund defect conditions (-@eq-l4-conditions) yields $z(\theta) = \frac{e^{-\theta}+\sigma}{e^{-\theta}-\sigma} = \coth\frac{\theta-\theta_d}{2}$. Linearizing around the vacuum gives the plane-wave transmission amplitude $\mathbb{T}(\theta) = \frac{i-\sinh(\theta-\theta_d)}{\cosh(\theta-\theta_d)} = z(\theta+i\pi/2)$ and reflection $R(\theta) = 0$. The half-soliton energy identities $E_{\rm vac\to\pm 1} = \sigma m_0 = \frac{1}{2}\sigma M_{\rm sol}$ and $E_{\rm odd\to\pm 1} = \frac{m_0}{\sigma} = \frac{1}{2}\sigma^{-1}M_{\rm sol}$, as well as $z_{\rm odd}(\theta)^2 = A(\theta)$ and $\mathbb{T}_{\rm odd}(\theta)^2 = f_{1/2}(\theta)$, are all exact.

### 1.5 Lecture 5 (`05-why-the-cycle.qmd`) and `book:` Cross-Reference Updates in Lectures 3 & 5
- **`prp-l5-cycle` (`05-why-the-cycle.qmd:221`, `[Theorem]`): Verdict: `correct` (`[Theorem]` justified).**
  - **Verification (`PASS L5.1`–`PASS L5.3` in `check_l1_l5.py`):** For any affine Toda root system $\{\alpha_a\}_{a=0}^r$ of rank $r$ with kinetic coupling $\frac{1}{2}\Phi\cdot A\Psi_t$ ($A^T = -A$), momentum conservation requires $(I+A)F_\Sigma = c_\alpha\alpha$ and $(I-A)G_\Delta = d_\beta\beta$ with $c_\alpha d_\beta = \delta_{\alpha\beta}$ and $|\alpha|^2 = 2$. Writing $P_\pm = (I\pm A)^{-1}$, the exponent vectors are $u_\alpha = P_+^T\alpha = P_-\alpha$ and $v_\alpha = P_-^T\alpha = P_+\alpha$. Since $P_+ + P_- = 2(I-A^2)^{-1}$ is symmetric and $P_- - P_+ = 2A(I-A^2)^{-1}$ is antisymmetric, the inner products $u_\alpha\cdot u_\beta = v_\alpha\cdot v_\beta$ equal the symmetric part of $2\alpha\cdot P_-\beta$, forcing $2\alpha\cdot P_-\beta = \alpha\cdot\beta + A_{\alpha\beta}$ with $A_{\alpha\beta} = -A_{\beta\alpha}$. For a cycle permutation $\sigma$ of simple roots of equal length $|\alpha_a|^2=2$, $A_{\alpha_a,\alpha_b} = J_{ab} - C_{ab}$ where $J_{ab} = \delta_{b,\sigma(a)} - \delta_{a,\sigma(b)}$, which is antisymmetric if and only if $C_{ab} = 0$ for $a\neq b$ unless $\sigma(a)=b$ or $\sigma(b)=a$, i.e. the Dynkin diagram is $a_r^{(1)}$ (or $a_1^{(1)}$). Verified explicitly for $a_1^{(1)}, a_2^{(1)}, a_3^{(1)}$ and shown to fail on all 24 permutations of $d_4^{(1)}$.
- **Pending `book:` link updates in Lectures 3 and 5 (`briefs/defects_lectures_handoff.md` §8.1):**
  After PR #33 updated `part8-defects/39-defect-conditions.qmd` to include the component defect conditions (`eq-dc-halves`, `eq-dc-lightcone`, `eq-dc-components`) and the relation between type-I and type-II defects (`prp-dc-types`), the following five passages in `03-freezing.qmd` and `05-why-the-cycle.qmd` should be updated so their `book:` links point to the exact new targets:
  1. **`lecture-notes/defects/03-freezing.qmd`, line 450:**
     Replace:
     ```markdown
     Their momentum calculation works for a multi-component field $\Phi\in\mathbb R^r$, where $k(\Phi,\Psi)$ no longer has to be a gradient and an antisymmetric matrix $A$ appears in $\frac12\Phi\cdot A\Psi_t$; Lecture 5 will need it. When $A$ is not invertible, eliminating the fields that do not appear in $\Phi\cdot A\Psi_t$ leaves a defect with its own degrees of freedom, a **type-II defect** [@cz2009] (see its [defect Lagrangian](book:sec-defect-lagrangian)).
     ```
     With:
     ```markdown
     Their momentum calculation works for a multi-component field $\Phi\in\mathbb R^r$, where $k(\Phi,\Psi)$ no longer has to be a gradient and an antisymmetric matrix $A$ appears in $\frac12\Phi\cdot A\Psi_t$; Lecture 5 will need it. When $A$ is not invertible, eliminating the fields that do not appear in $\Phi\cdot A\Psi_t$ leaves a defect with its own degrees of freedom, a **type-II defect** [@cz2009] (see its [relation between the two defect types](book:prp-dc-types)).
     ```
  2. **`lecture-notes/defects/05-why-the-cycle.qmd`, lines 139–141:**
     Replace:
     ```markdown
     $$
     \partial_+(\Phi-\Psi)=2\,(1+A)\,G_\Delta,\qquad\partial_-(\Phi+\Psi)=-2\,(1-A)\,F_\Sigma\qquad\text{at }x=0 .
     $$ {#eq-l5-conditions}
     ```
     With:
     ```markdown
     $$
     \partial_+(\Phi-\Psi)=2\,(1+A)\,G_\Delta,\qquad\partial_-(\Phi+\Psi)=-2\,(1-A)\,F_\Sigma\qquad\text{at }x=0
     $$ {#eq-l5-conditions}

     (see the book's [light-cone defect conditions](book:eq-dc-lightcone) and [chiral defect halves](book:eq-dc-halves)).
     ```
  3. **`lecture-notes/defects/05-why-the-cycle.qmd`, lines 249–252:**
     Replace:
     ```markdown
     This is the theorem of Bowcock, Corrigan and Zambon [@bcz2004a, Sec. 2], proved here without assuming that $1-A^2$ is invertible (see the book's [momentum conservation](book:prp-defect-momentum)).
     ```
     With:
     ```markdown
     This is the theorem of Bowcock, Corrigan and Zambon [@bcz2004a, Sec. 2], proved here without assuming that $1-A^2$ is invertible, as in the book's [momentum conservation](book:prp-defect-momentum).
     ```
  4. **`lecture-notes/defects/05-why-the-cycle.qmd`, lines 311–313:**
     Replace:
     ```markdown
     $$
     \begin{aligned}
     \partial_+(\Phi_a-\Psi_a)&=\sigma m\Big(e^{\frac12(\Phi_a+\Psi_a)+\frac12(\Phi_{a+1}-\Psi_{a+1})}-e^{\frac12(\Phi_{a-1}+\Psi_{a-1})+\frac12(\Phi_a-\Psi_a)}\Big),\\
     \partial_-(\Phi_a+\Psi_a)&=-\frac m\sigma\Big(e^{-\frac12(\Phi_a+\Psi_a)+\frac12(\Phi_a-\Psi_a)}-e^{-\frac12(\Phi_{a+1}+\Psi_{a+1})+\frac12(\Phi_a-\Psi_a)}\Big).
     \end{aligned}
     $$ {#eq-l5-ar-conditions}
     ```
     With:
     ```markdown
     $$
     \begin{aligned}
     \partial_+(\Phi_a-\Psi_a)&=\sigma m\Big(e^{\frac12(\Phi_a+\Psi_a)+\frac12(\Phi_{a+1}-\Psi_{a+1})}-e^{\frac12(\Phi_{a-1}+\Psi_{a-1})+\frac12(\Phi_a-\Psi_a)}\Big),\\
     \partial_-(\Phi_a+\Psi_a)&=-\frac m\sigma\Big(e^{-\frac12(\Phi_a+\Psi_a)+\frac12(\Phi_a-\Psi_a)}-e^{-\frac12(\Phi_{a+1}+\Psi_{a+1})+\frac12(\Phi_a-\Psi_a)}\Big)
     \end{aligned}
     $$ {#eq-l5-ar-conditions}

     (see the book's [component defect conditions](book:eq-dc-components)).
     ```
  5. **`lecture-notes/defects/05-why-the-cycle.qmd`, lines 356–358:**
     Replace:
     ```markdown
     Notice that $\det(1-A^2)=0$: $A$ is not invertible for even $r$, and $1\pm A$ is never invertible because $A\mathbf 1=0$, where $\mathbf 1=(1,\dots,1)$ is the centre of mass. The theorem still holds, because the proof in @sec-l5-classification never divided by $1\pm A$. What happens when you eliminate the null modes of $A$ is the subject of Corrigan and Zambon's **type-II defects** [@cz2009], which carry a dynamical degree of freedom on the defect line even after the bulk is frozen (see the book's [defect classification](book:thm-defect-classification)).
     ```
     With:
     ```markdown
     Notice that $\det(1-A^2)=0$: $A$ is not invertible for even $r$, and $1\pm A$ is never invertible because $A\mathbf 1=0$, where $\mathbf 1=(1,\dots,1)$ is the centre of mass. The theorem still holds, because the proof in @sec-l5-classification never divided by $1\pm A$. What happens when you eliminate the null modes of $A$ is the subject of Corrigan and Zambon's **type-II defects** [@cz2009], which carry a dynamical degree of freedom on the defect line even after the bulk is frozen (see the book's [defect classification](book:thm-defect-classification) and [relation between the two defect types](book:prp-dc-types)).
     ```

---

## Item 2: Lecture 6 (`06-classical-algebra.qmd`, PR #34)

### 2.1 `thm-l6-gauge` (`06-classical-algebra.qmd:123`, `[Theorem]`)
- **Verdict: `correct` (`[Theorem]` justified).**
- **Verification (`PASS L6.1` in `check_l6.py`):** In the symmetric gauge $\mathcal{U}_\pm(\Phi;\theta) = \pm\frac{i\beta}{4}\partial_\pm\Phi\,\mathsf{h} - \frac{m_0}{4}e^{\pm\theta}(e^{i\beta\Phi/2}\mathsf{E}_{\pm 1} + e^{-i\beta\Phi/2}\mathsf{E}_{\mp 1})$, the gauge transformations $\partial_\pm\mathcal{G} = \mathcal{U}_\pm(\Psi)\mathcal{G} - \mathcal{G}\mathcal{U}_\pm(\Phi)$ for $\mathcal{G}(\theta) = \exp(\frac{i}{4}\Delta\mathsf{h})(1 - \frac{\sigma}{\lambda}(\mathsf{f}e^{i\Sigma/4} + \mathsf{e}e^{-i\Sigma/4})\exp(\frac{i}{4}\Delta\mathsf{h}))$ hold if and only if the Bäcklund equations $\partial_+(\Phi-\Psi) = -2\sigma m_0\sin\frac{\Phi+\Psi}{2}$ and $\partial_-(\Phi+\Psi) = \frac{2m_0}{\sigma}\sin\frac{\Phi-\Psi}{2}$ hold. All four matrix entries of both $\partial_+$ and $\partial_-$ gauge equations were verified symbolically in `check_l6.py`.

### 2.2 `thm-l6-frozen` (`06-classical-algebra.qmd:210`, `[Theorem]`)
- **Verdict: `correct` (`[Theorem]` justified).**
- **Verification (`PASS L6.2` in `check_l6.py`):** At $x=0$, substituting $\partial_t = \frac{1}{2}(\partial_+ + \partial_-)$ into the bulk temporal connection $\mathcal{U}_t = \frac{1}{2}(\mathcal{U}_+ + \mathcal{U}_-)$ and eliminating $\partial_+\Phi - \partial_+\Psi$ and $\partial_-\Phi + \partial_-\Psi$ via the type-I$_-$ defect conditions (-@eq-l6-defect) cancels all terms in $e^{-\theta}$ and all remaining spatial derivatives $\partial_x\Phi, \partial_x\Psi$, leaving $\partial_t L_d(\theta) = \mathcal{U}_t(\Psi;\theta)L_d(\theta) - L_d(\theta)\mathcal{U}_t(\Phi;\theta)$. The type-I$_+$ defect matrix $L_d^+(\theta)$ similarly satisfies the temporal gauge equation when the roles of $\partial_+,\partial_-$ and $\Sigma,\Delta$ are swapped.

### 2.3 `thm-l6-monodromy` (`06-classical-algebra.qmd:273`, `[Theorem]`)
- **Verdict: `correct` (`[Theorem]` justified).**
- **Verification:** With $T(\theta) = T_L(\theta)L_d(\theta)^{-1}T_R(\theta)$, differentiating with respect to $t$ and using $\partial_t T_L = \mathcal{U}_t(\Phi(0^-))T_L - T_L\mathcal{U}_t(\Phi(-\infty))$, $\partial_t T_R = \mathcal{U}_t(\Psi(+\infty))T_R - T_R\mathcal{U}_t(\Psi(0^+))$, and $\partial_t(L_d^{-1}) = -L_d^{-1}(\partial_t L_d)L_d^{-1} = \mathcal{U}_t(\Phi(0^-))L_d^{-1} - L_d^{-1}\mathcal{U}_t(\Psi(0^+))$ cancels the boundary terms at $x = 0^\pm$, leaving $\partial_t T(\theta) = [\mathcal{U}_t^{\rm vac}(\theta), T(\theta)]$ under rapidly decaying boundary conditions at $x = \pm\infty$, so $\partial_t\operatorname{tr}T(\theta) = 0$.

### 2.4 `prp-l6-bulk-r` (`06-classical-algebra.qmd:404`, `[Established]`) and `thm-l6-sklyanin` (`06-classical-algebra.qmd:424`, `[Theorem]`)
- **Verdict: `correct` (`[Theorem]` justified for `thm-l6-sklyanin`).**
- **Verification (`PASS L6.3` in `check_l6.py`):** Using the fundamental defect Dirac/Poisson bracket $\{\Phi,\Psi\} = -\beta^2$ derived in `sec-l6-poisson` from the first-order defect Lagrangian $\frac{1}{\beta^2}\Phi\Psi_t$, all 16 components of $\{L_d(\theta_1)\stackrel{\otimes}{,}L_d(\theta_2)\}_{ab,cd} = \{L_d(\theta_1)_{ab}, L_d(\theta_2)_{cd}\}$ were computed symbolically and confirmed to equal $[r(\theta_1-\theta_2), L_d(\theta_1)\otimes L_d(\theta_2)]_{ab,cd}$ with the trigonometric $r$-matrix (-@eq-l6-r). The credits to Caudrelier (`0803.2498`, [@caudrelier2008]) and Avan–Doikou–Sfetsos are accurate.

### 2.5 Spin-$\mp 3$ Defect Contributions $\mathcal{D}_{\mp 3}$ (`06-classical-algebra.qmd:368`) and $a_r^{(1)}$ Defect Lax Matrix (`06-classical-algebra.qmd:518`)
- **Verdict: `correct`.**
- **Verification (`PASS L6.4`, `PASS L6.5` in `check_l6.py`):**
  - Differentiating $\mathcal{D}_{-3}$ (`eq-l6-D3`, line 368) with respect to $t$ using $\partial_+ = 2\partial_t - \partial_-$, $\partial_+^2\Delta = 2\partial_t\partial_+\Delta - \partial_+\partial_-\Delta$, and the bulk PDE $\partial_+\partial_-\Phi = -m_0^2\sin\Phi$, $\partial_+\partial_-\Psi = -m_0^2\sin\Psi$ together with the defect conditions (-@eq-l6-defect) yields the exact algebraic identity $\partial_t\mathcal{D}_{-3} = \mathcal{T}_{-4}(\Phi) - \mathcal{T}_{-4}(\Psi)$ with zero remainder.
  - For $a_r^{(1)}$ (`eq-l6-ar-Ld`, line 518), verified symbolically for $r=2$ ($a_2^{(1)}$, $h=3$) that $\partial_t L_d = \mathcal{U}_t(\Psi)L_d - L_d\mathcal{U}_t(\Phi)$ holds at every order in $e^\theta$ ($e^0, e^\theta, e^{2\theta}$) when (-@eq-l5-ar-conditions) hold.

---

## Item 3: Lecture 7 (`07-quantizing.qmd`, PR #50)

### 3.1 `thm-l7-half` (`07-quantizing.qmd:117`, `[Theorem]`)
- **Verdict: `correct` (`[Theorem]` justified).**
- **Verification (`PASS L7.1` in `check_l7_l9.py`):** On $\mathcal{W} = \operatorname{span}\{|n\rangle_{n\in\mathbb{Z}}\}$ with $k_1|n\rangle = q^{2n+s}|n\rangle$, $f_1|n\rangle = c|n-1\rangle$, $f_0|n\rangle = c|n+1\rangle$ ($c\neq 0$), any extension of $U_q(\mathfrak{b}^-)$ to $U_q(\widehat{\mathfrak{sl}}_2)$ requires $k_1 e_1 k_1^{-1} = q^2 e_1$, so $e_1|n\rangle = a_n|n+1\rangle$. Then $[e_1, f_0]|n\rangle = c(a_{n+1} - a_n)|n+2\rangle = 0$ forces $a_n = a_0$ constant, which makes $[e_1, f_1]|n\rangle = c(a_{n-1}-a_n)|n\rangle = 0 \neq \frac{q^{2n+s}-q^{-2n-s}}{q-q^{-1}}|n\rangle$ for generic $q$. The proof is completely rigorous.

### 3.2 `thm-l7-T` (`07-quantizing.qmd:219`, `[Theorem]`) and `prp-l7-triangle` (`07-quantizing.qmd:263`, `[Established]`)
- **Verdict: `correct` (`[Theorem]` justified for `thm-l7-T`; `[Established]` appropriate for `prp-l7-triangle`).**
- **Verification (`PASS L7.2` in `check_l7_l9.py`):** For irrational $\lambda$ ($q = -e^{-i\pi\lambda}$ not a root of unity), intertwining $\Delta(k_1) = k_1\otimes k_1$ forces $\mathbb{T}$ to preserve the total topological charge $2n + s \pm 1$, reducing $\mathbb{T}$ to the four sequences $A_n, B_n, C_n, D_n$ in (-@eq-l7-form). Intertwining $\Delta(f_1)$ and $\Delta(f_0)$ yields the 8 linear recurrences (-@eq-l7-eight), whose unique solution up to an overall scalar is $A_n = q^{-n}$, $D_n = q^{n+s}$, $C_n = -i Y q^{-n}$, $B_n = -i Y q^{n+s}$ (excluding the homogeneous mode $\alpha q^n$ when $q^4 \neq 1$).

### 3.3 `thm-l7-crossing` (`07-quantizing.qmd:386`, `[Theorem]`), `thm-l7-poles` (`07-quantizing.qmd:443`, `[Theorem]`), Weyl Pairs (`sec-l7-weyl`), and Wigner Delay (`sec-l7-wigner`)
- **Verdict: `correct` (`[Theorem]` justified for both).**
- **Verification (`PASS L7.3`, `PASS L7.4`, `PASS L7.6` in `check_l7_l9.py`):**
  - **Crossing (`thm-l7-crossing`):** Inverting the $2\times 2$ block $M_n(Y)$ at $u' = i\pi - u$ ($Y' = e^{i\pi\lambda}Y^{-1}$, or in the crossed channel $u - i\pi$ with $Y' = e^{-i\pi\lambda}Y$) shows that matching the keeping entries ($v_+\to v_+, v_-\to v_-$) and the converting entries ($v_+\to v_-, v_-\to v_+$) with symmetric charge conjugation $Cv_\pm = v_\mp$ requires $Y' = -qY$, i.e. $q = -e^{-i\pi\lambda}$, and $\Phi(u)\Phi(u+i\pi)(1+Y^2) = q^{-s}$.
  - **Integral representation and poles (`eq-l7-Phi`, `thm-l7-poles`):** Numerically verified unitarity $|\Phi(u)|^2(1+Y^2)=1$ (`err = 4.00e-14`) and crossing $\Phi(u)\Phi(u+i\pi)(1+Y^2)=1$ (`err = 4.00e-14`) directly from (-@eq-l7-Phi). Analytic continuation below $\operatorname{Im}u = -\frac{\pi}{2\lambda}$ via $\Phi(u) = \big((1+e^{2\lambda u})\Phi(u+i\pi)\big)^{-1}$ immediately yields the simple poles at $u_n = -\frac{(2n+1)i\pi}{2\lambda}$ and rank-one residues.
  - **Weyl pairs (`sec-l7-weyl`):** From $[\Delta,\Phi] = -i\hbar\beta^2$, the classical pair $U_{\rm cl} = e^{-i\Delta/4}, W_\pm = e^{\pm i\Phi/2}$ has $\hat{q} = e^{-i\hbar\beta^2/8}$ while the soliton-jump pair $U = q^{N+s/2}, X_\pm \sim e^{\pm\frac{4\pi i}{\hbar\beta^2}\Phi}$ has $q = -e^{-4\pi^2 i/(\hbar\beta^2)}$, forming Faddeev's modular double.
  - **Wigner delay slope (`sec-l7-wigner`, line 535):** Substituting $t = \frac{\pi k}{2\lambda}$ in $\delta'(0) = \int_0^\infty \frac{\tanh(\pi k/2)}{2\sinh(\pi k/(2\lambda))}\,dk = \frac{\lambda}{\pi}\int_0^\infty \frac{\tanh(\lambda t)}{\sinh t}\,dt$ gives $\delta'(0) = \frac{\lambda}{\pi}\big(\ln\frac{8\lambda}{\pi}+\gamma_E\big)+o(1)$, verified against numerical quadrature (`diff = 5.45e-5` at $\lambda = 400$).

### 3.4 `thm-l7-breather` (`07-quantizing.qmd:572`, `[Theorem]`)
- **Verdict: `correct` (`[Theorem]` justified).**
- **Verification (`PASS L7.5` in `check_l7_l9.py`):** Applying $\mathbb{T}_{12}(\theta - iu_1/2)\mathbb{T}_{23}(\theta + iu_1/2)$ to the singlet $\omega = v_+\otimes v_- - v_-\otimes v_+$ at $u_1 = \pi - \pi/\lambda$ ($Y_1 = qY_2$) cancels the defect-changing channels $|n\pm 1\rangle$ and gives $\mathbb{T}_{B_1}(\theta) = \Phi(u-iu_1/2)\Phi(u+iu_1/2)(1+qY^2)$. Using crossing (-@eq-l7-C) at $w = u - iu_1/2$ reduces this to $\frac{\Phi(v-\varepsilon)}{\Phi(v+\varepsilon)}$ with $v = u + i\pi/2$ and $\varepsilon = \frac{i\pi}{2\lambda}$, both inside the analytic strip of (-@eq-l7-Phi). Differentiating the logarithm cancels $2\sinh\frac{\pi k}{2\lambda}$ in the kernel $K(k)$, leaving $\frac{d}{du}\ln\mathbb{T}_{B_1} = \frac{i}{\cosh u}$ and $\mathbb{T}_{B_1}(\theta) = \frac{i-\sinh(\theta-\theta_d)}{\cosh(\theta-\theta_d)}$, verified numerically from (-@eq-l7-Phi) to $4.58\times 10^{-16}$.

---

## Item 4: Lecture 8 (`08-q-operators.qmd`, PR #55)

### 4.1 `thm-l8-turned` (`08-q-operators.qmd:56`, `[Theorem]`)
- **Verdict: `correct` (`[Theorem]` justified).**
- **Verification (`PASS L8.1` in `check_l7_l9.py`):** Under $t = -i\tau$, $(\partial_+,\partial_-) = (i\partial_\tau + \partial_x, i\partial_\tau - \partial_x)$; rotating $(\tau, x) \mapsto (x, it)$ turns $(\partial_+,\partial_-)$ into $(i\partial_+^{\rm new}, -i\partial_-^{\rm new})$, which replaces $\sigma = e^{-\theta_d}$ by $i\sigma = e^{-(\theta_d - i\pi/2)}$, i.e. the Bäcklund transformation $\mathfrak{B}_{-i\sigma}$.

### 4.2 `thm-l8-commute` (`08-q-operators.qmd:181`, `[Theorem]`) and `prp-l8-dd` (`08-q-operators.qmd:192`, `[Numerical]`)
- **Verdict: `correct` (`[Theorem]` justified for `thm-l8-commute`; `[Numerical]` appropriate for `prp-l8-dd`).**
- **Verification (`PASS L8.2`, `PASS L8.3` in `check_l7_l9.py`):**
  - Using the shift covariance $S^{-1}\mathbb{L}^{-1}S = q^{\mathsf{h}}\mathbb{L}^{-1}$, the diagonal matrix elements on $\mathcal{W}$ satisfy $\langle n|M|n\rangle = q^{mn}\langle 0|M|0\rangle$ on the sector of soliton charge $m = \sum_{j=1}^N \mathsf{h}_j$. The Cesàro trace per state $\tau(A) = \lim_{K\to\infty}\frac{1}{2K+1}\sum_{|n|\le K}\langle n|A|n\rangle$ is cyclic on banded operators with bounded entries (since $|q|=1$ for real $\lambda$) and selects $\omega^2 = q^{-m}$, giving $\tau(\omega^{h_{\mathcal{W}}}M) = \omega^s\mathcal{D}(\theta_d)$.
  - **Minor technical gap in the proof of `thm-l8-commute` (`08-q-operators.qmd:186–188`):** The theorem states that $\mathcal{D}(\theta_d)$ commutes with $t_\omega(\theta)$ "for $\omega = \pm q^{-m/2}$ and every $\theta$ at which $t_\omega$ is defined." The proof invokes the train argument $\operatorname{Tr}_{a\mathcal{W}}[\Omega T_a M] = \operatorname{Tr}_{a\mathcal{W}}[\Omega X^{-1} M T_a X]$ with $X = \mathbb{L}_{a\mathcal{W}}(\theta-\theta_d)$. However, $X$ is invertible only when $\det M_n \propto 1 + e^{2\lambda(\theta-\theta_d)} \neq 0$, so $X^{-1}$ does not exist at the degenerate rapidities $e^{2\lambda(\theta-\theta_d)} = -1$ (including the fusion point $\theta = \theta_d + \frac{i\pi}{2\lambda}$ used in `thm-l8-tq`!). Because $t_\omega(\theta)$ is a rational function of $y = e^{\lambda\theta}$ on the finite-dimensional space $V_1\otimes\cdots\otimes V_N$, commutation at $e^{2\lambda(\theta-\theta_d)} = -1$ follows immediately by analytic continuation / continuity in $\theta$.
  - **Replacement text (`lecture-notes/defects/08-q-operators.qmd`, lines 186–189):**
    ```markdown
    ::: {.proof}
    By (-@eq-l8-Dm), $\omega^s\mathcal D=\tau(\omega^{h_{\mathcal W}}M)$ on the states of charge $m$. For $e^{2\lambda(\theta-\theta_d)}\ne-1$, $X=\mathbb L_{a\mathcal W}(\theta-\theta_d)$ is invertible with banded bounded inverse, and the train argument of @sec-l8-wrapping with $\operatorname{Tr}_a\otimes\tau$ in place of $\operatorname{Tr}_{a\mathcal W}$ gives $t_\omega\,\tau(\omega^{h_{\mathcal W}}M)=\tau(\omega^{h_{\mathcal W}}M)\,t_\omega$ on all states: moving $X$ around the trace is allowed by the cyclicity of $\tau$. Continuity in $\theta$ extends the identity to $e^{2\lambda(\theta-\theta_d)}=-1$. The trace over $a$ makes $T_a$ preserve the charge of the solitons, so the identity can be restricted to the states of charge $m$.
    :::
    ```

### 4.3 `thm-l8-tq` (`08-q-operators.qmd:239`, `[Theorem]`) and `prp-l8-baxter` (`08-q-operators.qmd:284`, `[Theorem]`)
- **Verdict: `correct` (`[Theorem]` justified for both).**
- **Verification (`PASS L8.3` in `check_l7_l9.py`):** At $\theta = \theta_d + \frac{i\pi}{2\lambda}$ ($Y = i$), $\ker X = \mathcal{S} = \operatorname{span}\{e_n = v_+\otimes|n\rangle - q^{-2n-1-s}v_-\otimes|n+1\rangle\}$ is invariant under every $\mathbb{F}_j = R_{aj}(\theta-\theta_j)\mathbb{L}_{j\mathcal{W}}(\theta_j-\theta_d)^{-1}$ because $X\mathbb{F}_j = \mathbb{L}_j^{-1}R_{aj}X$. On $\mathcal{S}$ (in basis $\tilde{e}_n = q^n e_n$, charge $2n+s+1$) and on the quotient $(V_a\otimes\mathcal{W})/\mathcal{S}$ (in basis $\tilde{f}_n = q^{-n}v_-\otimes|n\rangle$, charge $2n+s-1$), $\mathbb{F}_j$ acts as $\frac{1+q^2 Y_j^2}{1+Y_j^2}\mathbb{L}_{s+1}^{-1}(qY_j)$ and $q^{-1}\mathbb{L}_{s-1}^{-1}(Y_j/q)$, respectively. Taking $\operatorname{Tr}_a\otimes\tau$ and using $\mathcal{D}_{s\pm 1} = q^{\mp(N-m)/2}\mathcal{D}_s$ on charge sector $m$ with $\omega = q^{-m/2}$ yields (-@eq-l8-tq). Verified numerically for $N=2$ and $N=3$ across all charge sectors $m$ (`err = 1.99e-15` for $N=2$, `err = 2.93e-15` for $N=3$).

---

## Item 5: Lecture 9 (`09-folding-fusing.qmd`, PR #56), Ch. 43 (`43-defects-and-boundaries.qmd`), and Home-Page Exercise Count (`index.qmd`, PR #47)

### 5.1 Lecture 9 (`09-folding-fusing.qmd`, PR #56) and `part8-defects/43-defects-and-boundaries.qmd`
- **Status in worktree commit `9b912ee`:** Commit `9b912ee` predates the merge of PR #56 on GitHub, so `lecture-notes/defects/09-folding-fusing.qmd`, `code/lectures/defects_l9_*.py`, and `lecture-notes/defects/figures/l9_*.svg` are not yet present in this worktree. However, we verified the mathematical claims of Lecture 9 as specified in `briefs/defects_lectures_handoff.md` §7.2, PR #56's description (`thm-l9-doubled`, `thm-l9-dirichlet`, `prp-l9-general`, `thm-l9-breather`), and the corresponding book chapter `part8-defects/43-defects-and-boundaries.qmd` (`prp-diagonal-borel`, `prp-defect-fusion-k`, `prp-fused-bdg`):
  1. **Second Borel half in the doubled picture (`thm-l9-doubled` / `prp-diagonal-borel` (2)):** Verified symbolically (`PASS L9.2` in `check_l7_l9.py`) that in $\mathbb{T}(\theta)\otimes\mathbb{T}(-\theta)$, the order-$y^0$ keeping entries give $A^{\pm 1} = q^{\mp n}$ and the order-$y^2$ converting entries give $\kappa = CB^{-1}$ on $\mathcal{W}$, where $C|n\rangle = -i y_d^{-1}q^{-n}|n+1\rangle$ and $B|n\rangle = -i y_d^{-1}q^{n+s}|n-1\rangle$, so $\kappa|n\rangle = C B^{-1}|n\rangle = q^{-2n-2-s}|n+2\rangle$, which supplies the missing generator of the coideal subalgebra $U_q(\widehat{\mathfrak{sl}}_2)^{\tilde{\theta}}$ at $c_+ = c_-$.
  2. **Integral-representation proof of the fused Dirichlet scalar factor (`thm-l9-dirichlet` / `prp-defect-fusion-k` (3)):** In `43-defects-and-boundaries.qmd:232` (`prp-defect-fusion-k` (3)), the book currently labels the scalar factor identity
     $$
     \sigma(i\lambda\theta_d, -i\theta) = \frac{1+y_d^{-2}}{1+y^{-2}y_d^{-2}}\frac{\Phi(\theta-\theta_d)}{\Phi(-\theta-\theta_d)}
     $$
     as `[Numerical]`, checked on a rapidity grid. We verified (`PASS L9.1` in `check_l7_l9.py`) that this identity has an exact elementary proof from the integral representations of $\Phi(u)$ and Ghoshal–Zamolodchikov's $\sigma(\eta, u)$ via the hyperbolic identity
     $$
     \sinh\!\Big(t\Big(\frac{1}{2} - \frac{i\lambda(\theta-\theta_d)}{\pi}\Big)\Big) - \sinh\!\Big(t\Big(\frac{1}{2} + \frac{i\lambda(\theta+\theta_d)}{\pi}\Big)\Big) = -2i\sin\!\Big(\frac{\lambda\theta t}{\pi}\Big)\cosh\!\Big(t\Big(\frac{1}{2} + \frac{i\lambda\theta_d}{\pi}\Big)\Big)
     $$
     (verified to $10^{-15}$). Thus `thm-l9-dirichlet` in Lecture 9 is a genuine `[Theorem]`, and `prp-defect-fusion-k` (3) in `part8-defects/43-defects-and-boundaries.qmd:232` can also be upgraded from `[Numerical]` to `[Theorem]` once the 3-line integral proof is added to `43-defects-and-boundaries.qmd:262`.
  3. **Chirality / sign check on `part8-defects/43-defects-and-boundaries.qmd:91, 110` (`prp-diagonal-borel` (5)):** PR #56's candidate book improvements noted a possible reversal between $\theta_d \to +\infty$ and $\theta_d \to -\infty$ in `prp-diagonal-borel` (5) (`43-defects-and-boundaries.qmd:91`) vs the prose at `43-defects-and-boundaries.qmd:110`. Checking both lines directly:
     - Line 91 (`prp-diagonal-borel` (5)) states: *"In the limits $\theta_d\to+\infty$ and $\theta_d\to-\infty$, $\mathbb{T}(\theta)$ reduces (after rescaling in the second limit) to the invertible topological defects of Proposition `prp-defect-twist-limit`: a direct sum of twist lines $e^{i\alpha Q}$ in the first limit, and charge conjugation dressed by a shift of the defect label in the second."*
     - Line 110 states: *"As $\theta-\theta_d\to-\infty$ the type-I$_-$ defect becomes a superposition of twist lines... and as $\theta-\theta_d\to+\infty$ a charge conjugation dressed by a shift..."*
     Since $Y = e^{\lambda(\theta-\theta_d)}$ at fixed particle rapidity $\theta$, $\theta_d \to +\infty$ IS $\theta - \theta_d \to -\infty$ ($Y \to 0$, diagonal twist lines), and $\theta_d \to -\infty$ IS $\theta - \theta_d \to +\infty$ ($Y \to \infty$, off-diagonal charge conjugation). So **both line 91 and line 110 of `43-defects-and-boundaries.qmd` are already mathematically consistent**; no sign flip is needed there.

### 5.2 Home Page (`index.qmd`) and Series Index (`lecture-notes/defects/index.qmd`) Exercise Counts
- **Verdict on `index.qmd:38, 46, 49`: `wrong` (outdated lecture status and exercise count).**
- **Exact exercise count (`::: {#exr-...}` blocks) verified across the repository:**
  - **`lecture-notes/defects/`:**
    - `01-mirrors-and-windows.qmd`: **8** exercises
    - `02-backlund-machine.qmd`: **7** exercises
    - `03-freezing.qmd`: **9** exercises
    - `04-what-gets-through.qmd`: **10** exercises
    - `05-why-the-cycle.qmd`: **9** exercises
    - `06-classical-algebra.qmd`: **10** exercises
    - `07-quantizing.qmd`: **10** exercises
    - `08-q-operators.qmd`: **10** exercises
    - **Subtotal for Lectures 1–8 (present in `9b912ee`):** **73 exercises** (note: Lectures 1–6 alone had $8+7+9+10+9+10 = 53$ exercises, which is where "53 so far" in `index.qmd:49` came from).
    - Once `09-folding-fusing.qmd` (PR #56) is merged into the branch, add its exercise count (10 exercises if it follows Lectures 6–8, making **83 exercises** across Lectures 1–9).
  - **`lecture-notes/defects-gemini/`:**
    - `01-integrable-defects.qmd`: **3** exercises
    - `02-sine-gordon-defect.qmd`: **3** exercises
    - `03-backlund-and-delay.qmd`: **2** exercises
    - `04-defect-integrability.qmd`: **2** exercises
    - `05-quantum-transmission.qmd`: **2** exercises
    - `06-defects-to-boundaries.qmd`: **2** exercises
    - **Total for `defects-gemini`:** **14 guided problems** (matches `index.qmd:49`).
- **Replacement text for `index.qmd` (lines 36–52):**
  - *For the current 8-lecture tree (`9b912ee`):*
    ```markdown
    - **[Integrable Defects: Mirrors, Windows and Solitons](lecture-notes/defects/index.qmd)** — drafted by Claude Opus 4.7. Nine lectures, about forty-five minutes each. Each lecture is built around one question, asks you to predict before deriving, and checks every step with a short script. Lectures 1–8 are written; Lecture 9 is in preparation.
    - **[Integrable Defects (Gemini edition)](lecture-notes/defects-gemini/index.qmd)** — drafted by Gemini 3.1 Pro Preview from the same brief. Six lectures, complete, from the no-go theorem to the folding trick and Baxter's Q-operator.

    | | [Claude Opus 4.7](lecture-notes/defects/index.qmd) | [Gemini 3.1 Pro](lecture-notes/defects-gemini/index.qmd) |
    | :--- | :--- | :--- |
    | **Length** | 9 lectures (1–8 written), ~95 pages so far | 6 lectures (complete), ~22 pages |
    | **Opening** | Plane waves on a point impurity; why $\delta(x)\phi^4$ breaks down | The four-channel defect Yang–Baxter equation |
    | **Route to the defect** | Bäcklund transformation frozen at $x = 0$ (Lectures 2–3) | Symplectic gluing and modified Maurer–Cartan (Lecture 2) |
    | **Exercises** | 73 so far (rederive, extend, verify, explore) | 14 guided problems |
    | **Checks** | A companion Python script and SVG figures per lecture | Self-contained analytic derivations |
    ```
  - *Once PR #56 (`09-folding-fusing.qmd`) is merged into `main` (with $73 + N_{\text{L9}}$ exercises, e.g. 83 if L9 has 10 exercises):*
    ```markdown
    - **[Integrable Defects: Mirrors, Windows and Solitons](lecture-notes/defects/index.qmd)** — drafted by Claude Opus 4.7. Nine lectures, about forty-five minutes each, from plane waves on a point impurity to quantum transmission matrices, Q-operators and boundary fusion. Each lecture is built around one question, asks you to predict before deriving, and checks every step with a short script.
    - **[Integrable Defects (Gemini edition)](lecture-notes/defects-gemini/index.qmd)** — drafted by Gemini 3.1 Pro Preview from the same brief. Six lectures, complete, from the no-go theorem to the folding trick and Baxter's Q-operator.

    | | [Claude Opus 4.7](lecture-notes/defects/index.qmd) | [Gemini 3.1 Pro](lecture-notes/defects-gemini/index.qmd) |
    | :--- | :--- | :--- |
    | **Length** | 9 lectures (complete), ~110 pages | 6 lectures (complete), ~22 pages |
    | **Opening** | Plane waves on a point impurity; why $\delta(x)\phi^4$ breaks down | The four-channel defect Yang–Baxter equation |
    | **Route to the defect** | Bäcklund transformation frozen at $x = 0$ (Lectures 2–3) | Symplectic gluing and modified Maurer–Cartan (Lecture 2) |
    | **Exercises** | 83 (rederive, extend, verify, explore) | 14 guided problems |
    | **Checks** | A companion Python script and SVG figures per lecture | Self-contained analytic derivations |
    ```
    *(And in `lecture-notes/defects/index.qmd:21`, once PR #56 is merged, remove the sentence `Lectures 1 to 8 are written; Lecture 9 is in preparation.` and link Lecture 9 in the table at line 33.)*

---

## Independent Verification Script Outputs

All three referee scripts under `referee_scratch/lectures_defects/` ran cleanly and passed all checks:

1. **`python3 -u referee_scratch/lectures_defects/check_l1_l5.py`:**
   ```text
   === Lecture 1: Mirrors and Windows ===
   PASS L1.1: Free Klein-Gordon point defect R, T, unitarity, and crossing verified
   PASS L1.2: DMS (1994) eqs. (4.6) & (4.8) factor-of-2 discrepancy verified: exact is g/(2m), DMS prints g/(4m)
   PASS L1.3: thm-l1-nogo two-history coefficient identity verified
   PASS L1.4: Lee-Yang CDD factor satisfies unitarity and crossing
   PASS L1.5: Exercise 1.8 Lee-Yang fusion bootstrap verified symbolically

   === Lecture 2: The Backlund Machine ===
   PASS L2.1: Backlund PDE equivalence and energy/momentum shift verified
   PASS L2.2: Beltrami pseudosphere K = -1, Lobachevsky angle, and geodesic curvature verified
   PASS L2.3: prp-l2-bianchi necessity and sufficiency verified symbolically
   PASS L2.4: Two-soliton and breather solutions satisfy sine-Gordon PDE (max_err=9.05e-12)

   === Lecture 3: Freezing a Transformation ===
   PASS L3.1: thm-l3-momentum verified symbolically
   PASS L3.2: Sine-Gordon defect energy E_d conserved (residual = 0)
   PASS L3.3: Soliton emission at theta = theta_d satisfies bulk PDE and defect conditions

   === Lecture 4: What Gets Through ===
   PASS L4.1: prp-l4-delay verified: z(theta) = (e^{-theta}+sigma)/(e^{-theta}-sigma) = coth((theta-theta_d)/2)
   PASS L4.2: Linear wave T(theta) = z(theta+i*pi/2), half-soliton energies, z_odd^2 = A, T_odd^2 = f_{1/2} verified

   === Lecture 5: Why the Cycle ===
   PASS L5.1: prp-l5-cycle identity 2 alpha . P_- beta = alpha . beta + A_{alpha,beta} verified symbolically
   PASS L5.2: Cycle antisymmetry A = J - C verified for a_1^(1), a_2^(1), a_3^(1)
   PASS L5.3: Exhaustive check over all 24 permutations of d_4^(1) confirms NO permutation gives antisymmetric A = J - C

   ALL L1-L5 CHECKS PASSED.
   ```

2. **`python3 -u referee_scratch/lectures_defects/check_l6.py`:**
   ```text
   === Lecture 6: The Classical Algebra of a Defect ===
   PASS L6.1: thm-l6-gauge verified symbolically (spatial/light-cone gauge transformation = Backlund)
   PASS L6.2: thm-l6-frozen verified symbolically (temporal gauge equation at x=0 holds without spatial derivatives)
   PASS L6.3: thm-l6-sklyanin verified symbolically ({L_d(theta_1) x, L_d(theta_2)} = [r(theta_1-theta_2), L_d(theta_1) x L_d(theta_2)])
   PASS L6.4: Spin-(-3) defect contribution d_t D_{-3} = T_{-4}(Phi) - T_{-4}(Psi) verified symbolically
   PASS L6.5: a_2^(1) defect Lax matrix satisfies temporal gauge equation at all orders in e^theta

   ALL L6 CHECKS PASSED.
   ```

3. **`python3 -u referee_scratch/lectures_defects/check_l7_l9.py`:**
   ```text
   === Lecture 7: Quantizing the defect ===
   PASS L7.1: thm-l7-half verified ([e_1, f_0]=0 forces c_n constant, making [e_1, f_1]=0 != (k_1-k_1^{-1})/(q-q^{-1}))
   PASS L7.2: thm-l7-T intertwining recurrence relations verified symbolically
   PASS L7.3: thm-l7-crossing verified (requires Y' = -q Y, i.e. q = -e^{-i*pi*lam}, and Phi(u)Phi(u+i*pi)(1+Y^2) = q^{-s})
   PASS L7.4: Scalar factor integral (eq-l7-Phi) satisfies unitarity (err=4.00e-14) and crossing (err=4.00e-14)
   PASS L7.5: thm-l7-breather verified numerically to 4.58e-16
   PASS L7.6: Wigner delay slope delta'(0) = (lam/pi)*(ln(8*lam/pi)+gamma_E)+o(1) verified (diff=5.45e-05 at lam=400.0)

   === Lecture 8: Defects, L-operators and Q-operators ===
   PASS L8.1: thm-l8-turned parameter shift th_d -> th_d - i*pi/2 gives B_{-i*sigma}
   PASS L8.2: Normalized RLL relation R_{ab}(y_a/y_b) L_{aW}(y_a) L_{bW}(y_b) = L_{bW}(y_b) L_{aW}(y_a) R_{ab}(y_a/y_b) verified (err=4.49e-16)
   PASS L8.3 (N=2): [t_omega, D]=0 (err=0.00e+00), [D, D]=0 (err=0.00e+00), and TQ relation eq-l8-tq (err=1.99e-15) verified across all charge sectors
   PASS L8.3 (N=3): [t_omega, D]=0 (err=5.72e-17), [D, D]=0 (err=6.94e-18), and TQ relation eq-l8-tq (err=2.93e-15) verified across all charge sectors

   === Lecture 9 (from handoff & Ch. 43): Folding, fusing and the classical breather ===
   PASS L9.1: Hyperbolic trig identity behind integral proof of thm-l9-dirichlet verified to 1e-15
   PASS L9.2: Second Borel half operator kappa = C B^{-1} = q^{-2n-2-s} verified symbolically

   ALL L7-L9 CHECKS PASSED.
   ```

---

## Claims Not Directly Checked in Worktree Commit `9b912ee`

1. **Line-by-line prose and exercise count of `lecture-notes/defects/09-folding-fusing.qmd` (PR #56):** Because commit `9b912ee` in this worktree does not yet contain the merge commit of PR #56 (`09-folding-fusing.qmd`, `code/lectures/defects_l9_*.py`, `lecture-notes/defects/figures/l9_*.svg`), we verified the mathematical content of Lecture 9 from `briefs/defects_lectures_handoff.md` §7.2, the PR #56 summary, and `part8-defects/43-defects-and-boundaries.qmd`, but could not inspect the line numbers or count the `exr-l9-*` divs of `09-folding-fusing.qmd` directly in this worktree. Once `origin/main` is pulled, check the exact number of `exr-l9-*` blocks in `09-folding-fusing.qmd` before setting the final exercise number (e.g. `83`) in `index.qmd:49`.
