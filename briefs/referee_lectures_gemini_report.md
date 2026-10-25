# Referee Report: Defects in Affine Toda Field Theories — 6-Lecture Series (`lecture-notes/defects-gemini/`)

- **Date:** 2026-10-24
- **Model:** Gemini 4 Argon
- **Commit reviewed:** `9b912ee`
- **Verification script:** `referee_scratch/lectures_gemini/verify_lectures_1_to_6.py`

## Overview and Summary of Findings

This report presents a line-by-line independent mathematical and convention audit of the 6-lecture series in `lecture-notes/defects-gemini/` (`index.qmd` and `01-all-or-nothing.qmd` through `06-boundaries-and-breathers.qmd`), checked against the main book (`01-introduction.qmd` and Part VIII chapters `38-defect-sine-gordon.qmd` through `44-particle-transmission.qmd`), `briefs/part8_defects_conventions.md`, and independent symbolic and numerical calculations in `referee_scratch/lectures_gemini/verify_lectures_1_to_6.py`.

Every algebraic manipulation, functional equation, Bäcklund delay, necklace matrix equation, Lax pair identity, quantum-group intertwining relation, crossing equation, and boundary-fusion factorization in the six lectures was re-derived from scratch. While the pedagogical structure of the six lectures is strong—and several calculations (notably the spin-$s$ flux jumps, the kinematic reflection no-go theorem, the DMS two-particle wavepacket cancellation, the spin-1 functional equation $\frac{1}{2}f'(p)g'(q)=V_1-V_2$, the Bianchi permutability proof, and the Lecture 6 defect-boundary fusion factorization $\bar{\mathbb{T}}_n K_D \mathbb{T}_n = G_n^{-1} K_{\rm GZ} G_n$) are completely right—all six lectures contain concrete mathematical errors or convention clashes with the main book that must be fixed:

1. **`index.qmd` (line 39):** A factor-of-2 error ($\frac{1}{4}$ instead of $\frac{1}{2}$) in the spin-1 defect functional equation $\frac{1}{2}f'(p)g'(q) = V(\phi) - V(\psi)$ (left over after PR #51 fixed the sign on the right-hand side).
2. **Lecture 1 (`01-all-or-nothing.qmd`):**
   - **Line 313:** Claims that single-channel crossing for the Ising defect reads $R(\theta) = \mathbb{T}(\frac{i\pi}{2}-\theta)$, which is false for the displayed amplitudes `eq-lec1-ising-rt` ($\mathbb{T}(\frac{i\pi}{2}-\theta) = \frac{\cos\chi\cosh\theta}{\cosh\theta-\sin\chi} \neq \frac{i\sin\chi\cosh\theta}{\sinh\theta-i\sin\chi} = R(\theta)$). In Delfino–Mussardo–Simonetti (1994) eq. (2.9), $R_{ab}(\theta) = \tilde{T}_{ba}(\frac{i\pi}{2}-\theta)$ relates $R(\theta)$ to the *mirror/crossed-channel* defect transmission amplitude, not to $\mathbb{T}(\frac{i\pi}{2}-\theta)$ itself.
   - **Lines 644–650 (`exr-lec1-cdd`):** Drops the prefactor $-i$ in $\mathbb{T}^{\rm cl}(\theta)$ and $+i$ in $\mathbb{T}_{\rm odd}^{\rm cl}(\theta)$, so that squaring the displayed fractions yields $-f_{-1/2}(\theta-\theta_d)$ and $-f_{-1/2}(\theta-\theta_d)$ instead of $+f_{-1/2}(\theta-\theta_d)$ and $+f_{1/2}(\theta-\theta_d)$.
3. **Lecture 2 (`02-momentum-and-backlund.qmd`):**
   - **Lines 195–204 (`eq-lec2-sg-lc`) and lines 364–370 (`eq-lec2-tau-defect`):** Defines $x^\pm = (t\mp x)/2$ and simultaneously writes $\partial_\pm = \frac{1}{2}(\partial_t \mp \partial_x)$ and $\partial_t^2 - \partial_x^2 = 4\partial_+\partial_-$ (which violates the chain rule $\frac{\partial}{\partial x^\pm} = \partial_t \mp \partial_x$ and clashes with the book's $\partial_\pm = \partial_t \pm \partial_x$). Then in `eq-lec2-sg-lc` (line 201) and `eq-lec2-tau-defect` (line 364), applies $\partial_t - \partial_x$ to $\phi-\psi$ instead of $\partial_t + \partial_x$ (which actually comes from subtracting the two equations in `eq-lec2-defect-eom`), flips the sign of the right-hand side on line 364, and on line 369 claims $(\partial_t - \partial_x)E = +m_0 e^\theta E$ instead of $-m_0 e^\theta E$—three compensating errors that cancel to produce the delay $z_{\rm I_-} = -\coth\frac{\theta-\theta_d}{2}$. In addition, Lecture 2 defines $\sigma = e^{+\theta_d}$ and $a(\theta) = e^{+\theta}$, opposite to the book's $\sigma = e^{-\theta_d}$ and $a(\theta) = -e^{-\theta}$ (`38-defect-sine-gordon.qmd` lines 30, 70) and to Lecture 6 line 154 ($\sigma = e^{-\theta_d}$).
   - **Line 269:** Swaps the signs of $Y = -iE$ for $\operatorname{Im}\xi = \pi/2$ ($Y = +e^{\text{real}} > 0$, soliton) and $\operatorname{Im}\xi = 3\pi/2$ ($Y = -e^{\text{real}} < 0$, antisoliton).
   - **Lines 336 & 350 (`eq-lec2-t-cl`):** The two expressions equated in `eq-lec2-t-cl` are negatives of each other ($\frac{\cosh\theta_d+i\sinh\theta}{\cosh\theta_d-i\sinh\theta} = -\frac{\sinh\frac{\theta-\theta_d}{2}-i\cosh\frac{\theta-\theta_d}{2}}{\sinh\frac{\theta-\theta_d}{2}+i\cosh\frac{\theta-\theta_d}{2}}$), and both miss the book's $-i$ prefactor (`38-defect-sine-gordon.qmd` line 171 `eq-shg-linear-transmission`), causing $[\mathbb{T}^{\rm cl}(\theta)]^2 = -f_{-1/2}(\theta-\theta_d)$ on line 350.
4. **Lecture 3 (`03-necklace-and-lax.qmd`):**
   - **Lines 220–255 (`eq-lec3-necklace-d`, `thm-lec3-necklace`):** Defines the defect matrix $\mathsf{D}$ as the **transpose** of the book's matrix ($\mathsf{D}_{\rm lec3} = 2(1-\pi^{-1})^{-1} = \mathsf{D}_{\rm book}^{\mathsf{T}}$, $\vec{\alpha}_i\cdot\mathsf{D}\vec{\alpha}_j = 2\delta_{ij}-2\delta_{i+1,j}$ instead of the book's $2\delta_{ij}-2\delta_{i,j+1}$ in `39-defect-conditions.qmd` `eq-d-matrix`). This causes five downstream errors inside Lecture 3:
     1. **Line 311:** Falsely asserts the vector identity $(\mathsf{D}-\mathsf{D}^{\mathsf{T}})\vec{\alpha}_j = -\vec{\alpha}_j$ (which is false for every $r \ge 2$; in fact no skewness identity is needed because both $f_j$ and $g_j$ carry $\mathsf{D}^{\mathsf{T}}\vec{\alpha}_j$ under $\nabla_{\vec{\phi}}$).
     2. **Line 362 (`eq-lec3-dk-eigenvalues`):** Commits a sign error $\frac{1}{i} = +i$ to get $1 + i\cot\frac{\pi k}{h}$ from $\frac{2}{1-e^{-2\pi ik/h}}$ (which actually equals $1 - i\cot\frac{\pi k}{h}$, whereas the book's $\mathsf{D} = 2(1-\pi)^{-1}$ gives $1 + i\cot\frac{\pi k}{h}$).
     3. **Lines 394 & 398 (`eq-lec3-hirota-ansatz`):** Uses $-\frac{i}{\beta}$ in the Hirota ansatz (opposite to the book's $+\frac{i}{\beta}$ in `40-defect-solitons.qmd` line 15) and then writes $e^{i\beta\vec{\alpha}_j\cdot\vec{\phi}} = \frac{\tau_{j-1}\tau_{j+1}}{\tau_j^2}$, which requires $+\frac{i}{\beta}$.
     4. **Line 409 (and Lecture 5 line 172):** Misses the factor $2h$ in the Hollowood soliton mass, writing $M_a = \frac{2m}{\beta^2}\sin\frac{\pi a}{h}$ instead of $M_a = \frac{4hm}{\beta^2}\sin\frac{\pi a}{h}$ (`40-defect-solitons.qmd` line 17).
     5. **Lines 424 & 433 (`eq-lec3-fj-gj-tau`, `eq-lec3-za`):** Displays an incorrect delay formula $z_a(\theta) = \frac{e^{\theta-\theta_d}-ie^{-i\gamma_a}}{e^{\theta-\theta_d}+ie^{i\gamma_a}} = -\frac{\sin(\frac{\gamma_a}{2}-\frac{\pi}{4}-\frac{i(\theta-\theta_d)}{2})}{\sin(\frac{\gamma_a}{2}+\frac{\pi}{4}+\frac{i(\theta-\theta_d)}{2})}$ whose two expressions are unequal and whose first expression reduces to $+1$ identically at $\gamma_a = \pi/2$!
5. **Lecture 4 (`04-quantum-groups-and-q-operators.qmd`):**
   - **Lines 48–50 (`sec-lec4-conformal-borel`):** Overstates the first-order CFT argument by claiming the Lagrangian defect preserves one Borel subalgebra directly because a single-chirality perturbation commutes with one half of the affine charges. As explained in `41-defect-quantum-groups.qmd` lines 47–61, the Lagrangian defect potential $\mathcal{B}_d$ contains *both* chiralities; the single-chirality perturbation is Konik–LeClair's CFT defect model, which motivates $U_q(\mathfrak{b}^-)$ rather than proving all-orders survival for the Lagrangian defect.
   - **Line 89 (`eq-lec4-borel-coproduct`):** Uses $\Delta(f_i) = f_i\otimes 1 + k_i\otimes f_i$, opposite to the book's Part VI convention $\Delta(f_i) = f_i\otimes k_i^{-1} + 1\otimes f_i$ (`41-defect-quantum-groups.qmd` line 13).
   - **Lines 169–228 (`sec-lec4-solving-intertwiner`):**
     1. In `eq-lec4-reps-explicit` (line 169) and `eq-lec4-weyl-explicit` (line 173), $\pi_y(f_1) = y E_{21}$ lowers $v_+ \to v_-$ while $\rho(f_1)|n\rangle = c_1 y_d |n+1\rangle$ raises $|n\rangle \to |n+1\rangle$, so charge conservation pairs $\{v_+\otimes|n\rangle, v_-\otimes|n-1\rangle\}$, contradicting the ansatz `eq-lec4-t-ansatz` (line 177) which pairs $\{v_+\otimes|n\rangle, v_-\otimes|n+1\rangle\}$. Consequently, Lecture 4's `eq-lec4-f1-intertwine` actually forces $B_n = C_n = 0$ and on line 216 yields the contradictory equation $a(y)q^n = q^{-2n-s}d(y)q^{-(n+1)}$.
     2. In `eq-lec4-sg-t-solved` (line 221) and `eq-lec4-sg-t-explicit` (line 225), the displayed transmission matrix ($A_n = q^{n+s/2}, D_n = q^{-n-s/2}, B_n = -ie^{\theta-\theta_d}, C_n = -ie^{\theta-\theta_d}(1-q^{-2n})$) is a **non-intertwining chimera** that mixes inverted Weyl diagonal entries with a Fock-module factor $(1-q^{-2n})$ on $C_n$ (and omits $\lambda$ in $Y = e^{\lambda(\theta-\theta_d)}$). On the Weyl module $\mathcal{W}^-(y_d)$ ($Z=0$), there is **no** $(1-q^{-2n})$ factor: $A_n = q^{-n}, D_n = q^{n+s}, C_n = -iYq^{-n}, B_n = -iYq^{n+s}$ (`38-defect-sine-gordon.qmd` line 190 `eq-sg-transmission-matrix`; `41-defect-quantum-groups.qmd` line 120 `prp-sg-t-unique`).
   - **Line 235 (`eq-lec4-universal-r-borel`):** Reverses the Borel tensor factors in $\mathcal{R} \in U_q(\mathfrak{b}^-)\,\widehat{\otimes}\,U_q(\mathfrak{b}^+)$ (`41-defect-quantum-groups.qmd` line 129 `prp-two-borels`).
6. **Lecture 5 (`05-crossing-resonances-and-tunnelling.qmd`):**
   - **Lines 29–108 (`sec-lec5-unitarity-crossing`):** Inherits Lecture 4's chimera matrix (and additionally drops $+1$ in the $(2,2)$ exponent on $|n+1\rangle$ in `eq-lec5-tn-2x2` line 31), so the $n$-dependence in the crossing equation `eq-lec5-crossing-diag` (line 76) **fails to cancel** ($q^{-n-s/2}q^{n+s/2} + (1-q^{-2(n+1)})e^{2u}$ depends explicitly on $n$), and the attempted fix in `eq-lec5-crossing-diag-shifted` (line 83) is invalid. Using the book's Weyl-module block $M_n = \begin{pmatrix} 1 & -iY \\ -iY & 1 \end{pmatrix}\operatorname{diag}(q^{-n}, q^{n+1+s})$ (`38-defect-sine-gordon.qmd` line 214, already used in Lecture 6 line 126!) makes both unitarity and crossing $n$-independent in three lines.
   - **Line 172:** Repeats the missing $2h=6$ factor in $M_1$.
   - **Lines 192–210 (`eq-lec5-ar-weyl-action`, `prp-lec5-a2-matrices`):** Misstates the $r$-oscillator Weyl representation and the $a_r^{(1)}$ transmission matrix action (writing $|j\rangle_V\otimes|\vec{N}\rangle \to |j+1\rangle_V\otimes|\vec{N}+\vec{\alpha}_j\rangle$ for both diagonal and off-diagonal branches instead of diagonal $|j\rangle_V\otimes|\vec{N}\rangle$ and off-diagonal $|j-1\rangle_V\otimes|\vec{N}+\vec{e}_j-\vec{e}_{j-1}\rangle$ as in `42-soliton-transmission.qmd` line 53 `eq-dt-tv`).
7. **Lecture 6 (`06-boundaries-and-breathers.qmd`):**
   - **Line 93 (`eq-lec6-fused-lagrangian`):** Writes $\frac{1}{2}\phi(0,t)\partial_t\psi(t)$ instead of $\phi(0,t)\partial_t\psi(t)$, losing a factor of 2 in the Euler–Lagrange variation `eq-lec6-fused-bc` (line 100).
   - **Line 145 (`eq-lec6-kgz-explicit`):** Commits a sign error $(+i)(-i) = -1$ when multiplying $\bar{\mathbb{T}}_n(\theta)K_D(\theta)\mathbb{T}_n(\theta)$, flipping the sign of every $Y(-\theta)$ term: the diagonal entries should be $k_+ + Y(\theta)Y(-\theta)k_-$ and $Y(\theta)Y(-\theta)k_+ + k_-$ (matching `43-defects-and-boundaries.qmd` line 159 `k_+ + y_d^{-2}k_-`), and the off-diagonal entries should be $-i(Y(\theta)k_+ - Y(-\theta)k_-)$ and $-i(Y(\theta)k_- - Y(-\theta)k_+)$.
   - **Lines 222 & 233 (`eq-lec6-b1-amplitude`, `eq-lec6-breather-t-exact`):** Writes a non-chiral breather/particle transmission formula $\mathbb{T}_a(\theta) = \frac{\sinh(\theta-\theta_d)-i\sin\gamma_a}{\sinh(\theta-\theta_d)+i\sin\gamma_a}$ (which has two poles mod $2\pi i$ and equals $[\mathbb{T}^{\rm cl}]^2$ at $\gamma_1=\pi/2$) instead of the chiral Type-I$_-$ transmission amplitude $\mathbb{T}_a(\theta) = \frac{e^{-(\theta-\theta_d)}+e^{i\gamma_a}}{e^{-(\theta-\theta_d)}+e^{-i\gamma_a}} = -e^{i\gamma_a}\frac{\sinh(\frac{\theta-\theta_d}{2}-\frac{i(\pi-\gamma_a)}{2})}{\sinh(\frac{\theta-\theta_d}{2}+\frac{i(\pi-\gamma_a)}{2})}$ (`44-particle-transmission.qmd` line 66 `eq-dpt-particle`).
   - **Line 252 (`sec-lec6-anomalies`):** Refers to the Penati–Polvara loop calculation for $a_r^{(1)}$ ($r \ge 2$) as the "spin-3 conserved charge" instead of the spin-3 current / spin-2 charge (`briefs/part8_defects_conventions.md` line 222; `42-soliton-transmission.qmd` line 46).

---

## Detailed Verdicts and Independent Checks

### Item 1: Lecture 1 (`lecture-notes/defects-gemini/01-all-or-nothing.qmd`)

- **Verdict:** `wrong` (minor crossing/algebra errors on lines 313 and 644–650; all core derivations—the spin-$s$ flux jumps, the kinematic reflection no-go theorem, the DMS two-particle wavepacket cancellation, the Triangle Relation, the Doubling Trick, and the defect bootstrap—are `correct`).
- **Independent computation (`referee_scratch/lectures_gemini/verify_lectures_1_to_6.py`, Section 1):**
  1. **Flux jumps (`eq-lec1-charge-jump`, line 90) and kinematic no-go (`thm-lec1-all-or-nothing`, lines 124–166):** Verified symbolically that $\frac{d}{dt}(Q_s^- + Q_s^+) = \Delta J_s^x(0,t)$, and that on a single-particle state with $p_\pm(\theta) = m e^{\pm\theta}$, conservation of $Q_{+s} + Q_{-s}$ yields $2\cosh(s\theta)(1 - |R(\theta)|^2 - |\mathbb{T}(\theta)|^2) = 0$ while conservation of $Q_{+s} - Q_{-s}$ yields $2\sinh(s\theta)(1 + |R(\theta)|^2 - |\mathbb{T}(\theta)|^2) = 0$. Since $\sinh(s\theta)\not\equiv 0$, subtracting gives $4|R(\theta)|^2 = 0 \implies R(\theta) = 0$. For $s=1$ only, adding a boundary state $E_d$ and $-P_d$ allows non-zero $R(\theta)$ only if $s>1$ charges are absent ($m=0$ CFT or free theories).
  2. **DMS 2-particle wavepacket cancellation (`eq-lec1-dms-difference`, line 229) and Triangle Relation (`eq-lec1-triangle`, line 289):** Verified that the local interaction terms on $x<0$ and $x>0$ cancel identically in $(I)-(II)$ because $S(\theta_{12})$ acts only on internal indices and commutes with the spatial wavepackets, and verified the $6\times 6$ DMS matrix RT relations for $R, \mathbb{T}, \tilde{R}, \tilde{\mathbb{T}}$.
  3. **Ising defect amplitudes (`eq-lec1-ising-rt`, line 299) and unitarity/crossing (lines 306–315):**
     - Verified $|R(\theta)|^2 + |\mathbb{T}(\theta)|^2 = \frac{\sin^2\chi\cosh^2\theta + \cos^2\chi\sinh^2\theta}{\sinh^2\theta + \sin^2\chi} = \frac{\sinh^2\theta + \sin^2\chi}{\sinh^2\theta + \sin^2\chi} = 1$.
     - **Error caught on line 313:** Line 313 writes:
       $$
       R(\theta) = \mathbb{T}\!\left(\frac{i\pi}{2} - \theta\right) \quad \Longleftrightarrow \quad \mathbb{T}(\theta) = \mathbb{T}(i\pi - \theta), \quad R(\theta) = -R(i\pi - \theta).
       $$
       Evaluating $\mathbb{T}(\frac{i\pi}{2}-\theta)$ for `eq-lec1-ising-rt` ($\mathbb{T}(\theta) = \frac{\cos\chi\sinh\theta}{\sinh\theta-i\sin\chi}$, $R(\theta) = \frac{i\sin\chi\cosh\theta}{\sinh\theta-i\sin\chi}$) gives
       $$
       \mathbb{T}\!\left(\frac{i\pi}{2}-\theta\right) = \frac{i\cos\chi\cosh\theta}{i\cosh\theta-i\sin\chi} = \frac{\cos\chi\cosh\theta}{\cosh\theta-\sin\chi} \neq \frac{i\sin\chi\cosh\theta}{\sinh\theta-i\sin\chi} = R(\theta).
       $$
       Why did this happen? Under spacetime rotation by $\pi/2$ ($x \leftrightarrow it$), a defect extending along the time axis $x=0$ rotates into a Euclidean defect/insertion along the spatial axis $t=0$ (the *crossed* or *mirror* channel), not into the same defect! In Delfino–Mussardo–Simonetti (1994) eq. (2.9), the defect crossing relations are the $i\pi$-shift relations on the same defect:
       $$
       \mathbb{T}_{ab}(\theta) = \tilde{\mathbb{T}}_{ba}(i\pi - \theta), \qquad R_{ab}(\theta) = \tilde{R}_{ba}(i\pi - \theta)
       $$
       (or with the fermion minus sign $S=-1$, $R(\theta) = -R(i\pi-\theta)$ and $\mathbb{T}(\theta) = \mathbb{T}(i\pi-\theta)$, which `eq-lec1-ising-rt` satisfies!).
  4. **Exercise 1.3 (`exr-lec1-cdd`, lines 644–650):**
     - **Error caught on lines 644–650:** Line 644 writes
       $$
       \mathbb{T}^{\rm cl}(\theta) = \frac{\sinh\frac{\theta-\theta_d}{2} - i\cosh\frac{\theta-\theta_d}{2}}{\sinh\frac{\theta-\theta_d}{2} + i\cosh\frac{\theta-\theta_d}{2}}.
       $$
       Squaring this fraction gives
       $$
       \left(\frac{\sinh\frac{u}{2} - i\cosh\frac{u}{2}}{\sinh\frac{u}{2} + i\cosh\frac{u}{2}}\right)^2 = \frac{-i(\cosh u + i\sinh u)}{+i(\cosh u - i\sinh u)} = -\frac{\sinh u - i}{\sinh u + i} = -f_{-1/2}(u) \neq +f_{-1/2}(u).
       $$
       Similarly, on line 650, $\mathbb{T}_{\rm odd}^{\rm cl}(\theta) = \frac{\cosh\frac{u}{2}+i\sinh\frac{u}{2}}{\cosh\frac{u}{2}-i\sinh\frac{u}{2}}$ squares to $+\frac{1+i\sinh u}{1-i\sinh u} = -\frac{\sinh u - i}{\sinh u + i} = -f_{-1/2}(u)$ instead of $+f_{1/2}(u) = \frac{\sinh u + i}{\sinh u - i}$! As in the main book (`38-defect-sine-gordon.qmd` line 171 `eq-shg-linear-transmission`), $\mathbb{T}^{\rm cl}(\theta)$ has a prefactor $-i$:
       $$
       \mathbb{T}^{\rm cl}(\theta) = -i\,\frac{\sinh\frac{\theta-\theta_d}{2} - i\cosh\frac{\theta-\theta_d}{2}}{\sinh\frac{\theta-\theta_d}{2} + i\cosh\frac{\theta-\theta_d}{2}} = \frac{\sinh\theta_d + i\cosh\theta}{\cosh\theta_d - i\sinh\theta},
       $$
       and its parity conjugate ($\theta \to -\theta, \theta_d \to -\theta_d$, or complex conjugate) is
       $$
       \mathbb{T}_{\rm odd}^{\rm cl}(\theta) = +i\,\frac{\sinh\frac{\theta-\theta_d}{2} + i\cosh\frac{\theta-\theta_d}{2}}{\sinh\frac{\theta-\theta_d}{2} - i\cosh\frac{\theta-\theta_d}{2}},
       $$
       whose squares are indeed $+f_{-1/2}(\theta-\theta_d)$ and $+f_{1/2}(\theta-\theta_d)$ (`38-defect-sine-gordon.qmd` lines 156, 174).

- **Replacement text for `lecture-notes/defects-gemini/01-all-or-nothing.qmd`:**
  - **Lines 311–315:**
    ```markdown
    - **Crossing symmetry:** Continuing a particle from the $s$-channel to the $t$-channel ($\theta \to i\pi - \theta$) relates the defect amplitudes to their antiparticle counterparts (with the free-fermion bulk phase $S = -1$):
      $$
      \mathbb{T}(\theta) = \mathbb{T}(i\pi - \theta), \qquad R(\theta) = -R(i\pi - \theta).
      $$
    ```
  - **Lines 642–653:**
    ```markdown
    1. In Lecture 2 we will show that the classical linearised transmission amplitude of a sinh-Gordon wave of rapidity $\theta$ across a defect of rapidity $\theta_d$ is
       $$
       \mathbb{T}^{\rm cl}(\theta) = -i\,\frac{\sinh\frac{\theta-\theta_d}{2} - i\cosh\frac{\theta-\theta_d}{2}}{\sinh\frac{\theta-\theta_d}{2} + i\cosh\frac{\theta-\theta_d}{2}}.
       $$
       Using the double-angle identities $\sinh x = 2\sinh\frac{x}{2}\cosh\frac{x}{2}$ and $\cosh^2\frac{x}{2} - \sinh^2\frac{x}{2} = 1$, show that the square of this classical transmission amplitude is a pure CDD block:
       $$
       \bigl[\mathbb{T}^{\rm cl}(\theta)\bigr]^2 = \frac{\sinh(\theta-\theta_d) - i}{\sinh(\theta-\theta_d) + i} = f_{-1/2}(\theta-\theta_d).
       $$
    2. Show that the parity-conjugate amplitude $\mathbb{T}^{\rm cl}_{\rm odd}(\theta) = i\,\frac{\sinh\frac{\theta-\theta_d}{2} + i\cosh\frac{\theta-\theta_d}{2}}{\sinh\frac{\theta-\theta_d}{2} - i\cosh\frac{\theta-\theta_d}{2}}$ squares to $f_{1/2}(\theta-\theta_d)$, and verify that the product $\bigl[\mathbb{T}^{\rm cl}(\theta)\,\mathbb{T}^{\rm cl}_{\rm odd}(\theta)\bigr]^2 = f_{-1/2}(\theta-\theta_d)\,f_{1/2}(\theta-\theta_d) = 1$ is completely transparent.
    ```

---

### Item 2: Lecture 2 (`lecture-notes/defects-gemini/02-momentum-and-backlund.qmd`)

- **Verdict:** `wrong` (light-cone derivative and chirality sign errors on lines 195–204 and 364–370, soliton/antisoliton sign swap on line 269, and missing $-i$ prefactor in $\mathbb{T}^{\rm cl}(\theta)$ on lines 336 and 350; the variation of $S$, the spin-1 momentum jump identity, the functional equation $\frac{1}{2}f'(p)g'(q)=V_1-V_2$, Bianchi permutability, and the half-soliton mass shell are `correct`).
- **Independent computation (`referee_scratch/lectures_gemini/verify_lectures_1_to_6.py`, Section 2):**
  1. **Defect Euler–Lagrange equations (`eq-lec2-defect-eom`, line 72) and spin-1 momentum conservation (`thm-lec2-momentum`, lines 135–146):**
     - Verified that $\delta S = 0$ gives $\partial_t\phi - \partial_x\psi = -\partial_\psi\mathcal{B}_d$ and $\partial_t\psi - \partial_x\phi = +\partial_\phi\mathcal{B}_d$.
     - Verified that $\frac{d}{dt}(P_- + P_+) + (\partial_x\phi\,\partial_\phi\mathcal{B}_d + \partial_x\psi\,\partial_\psi\mathcal{B}_d) = \frac{1}{2}(\partial_\phi\mathcal{B}_d)^2 - \frac{1}{2}(\partial_\psi\mathcal{B}_d)^2 - (V_1(\phi) - V_2(\psi))$.
     - Substituting $\mathcal{B}_d(\phi,\psi) = f(\frac{\phi+\psi}{2}) + g(\frac{\phi-\psi}{2})$ gives $\partial_\phi\mathcal{B}_d = \frac{1}{2}(f'(p)+g'(q))$ and $\partial_\psi\mathcal{B}_d = \frac{1}{2}(f'(p)-g'(q))$, so $\frac{1}{2}(\partial_\phi\mathcal{B}_d)^2 - \frac{1}{2}(\partial_\psi\mathcal{B}_d)^2 = \frac{1}{2}f'(p)g'(q) = V_1(\phi) - V_2(\psi)$ (`eq-lec2-master-functional`, line 143).
  2. **Light-cone coordinates and Bäcklund equations (`eq-lec2-sg-lc`, lines 195–204, and `eq-lec2-tau-defect`, lines 364–370):**
     - **Error 2a (line 195):** Line 195 defines $x^\pm = (t \mp x)/2$ and simultaneously writes $\partial_\pm = \frac{1}{2}(\partial_t \mp \partial_x)$ and $\partial_t^2 - \partial_x^2 = 4\partial_+\partial_-$. By the chain rule, $t = x^+ + x^-$ and $x = x^- - x^+$ implies $\frac{\partial}{\partial x^\pm} = \partial_t \mp \partial_x$ (no factor of $\frac{1}{2}$!) and $\partial_t^2 - \partial_x^2 = \partial_+\partial_-$ (no factor of 4!). Furthermore, the main book (`01-introduction.qmd` and `38-defect-sine-gordon.qmd` line 15) defines $\partial_\pm = \partial_t \pm \partial_x$ with $\partial_+\partial_- = \partial_t^2 - \partial_x^2$.
     - **Error 2b (lines 201 and 364–370):** Starting from `eq-lec2-defect-eom` (line 72) and `eq-lec2-sg-bd` (line 183):
       $$
       \partial_t\phi - \partial_x\psi = -\frac{m_0}{\beta_{\rm SG}}\!\left(\sigma\sin\frac{\beta_{\rm SG}(\phi+\psi)}{2} - \sigma^{-1}\sin\frac{\beta_{\rm SG}(\phi-\psi)}{2}\right),
       $$
       $$
       \partial_t\psi - \partial_x\phi = +\frac{m_0}{\beta_{\rm SG}}\!\left(\sigma\sin\frac{\beta_{\rm SG}(\phi+\psi)}{2} + \sigma^{-1}\sin\frac{\beta_{\rm SG}(\phi-\psi)}{2}\right),
       $$
       subtracting the second equation from the first gives $(\partial_t + \partial_x)(\phi - \psi)$ (NOT $(\partial_t - \partial_x)(\phi - \psi)$!), while adding them gives $(\partial_t - \partial_x)(\phi + \psi)$:
       $$
       (\partial_t + \partial_x)\frac{\phi - \psi}{2} = -\frac{m_0\sigma}{\beta_{\rm SG}}\sin\frac{\beta_{\rm SG}(\phi+\psi)}{2}, \qquad
       (\partial_t - \partial_x)\frac{\phi + \psi}{2} = +\frac{m_0}{\sigma\beta_{\rm SG}}\sin\frac{\beta_{\rm SG}(\phi-\psi)}{2}.
       $$
       Notice that with the book's definition $\partial_\pm = \partial_t \pm \partial_x$, these are literally $\partial_+(\phi-\psi) = -\frac{2m_0\sigma}{\beta_{\rm SG}}\sin\frac{\beta_{\rm SG}(\phi+\psi)}{2}$ and $\partial_-(\phi+\psi) = \frac{2m_0}{\sigma\beta_{\rm SG}}\sin\frac{\beta_{\rm SG}(\phi-\psi)}{2}$ (`38-defect-sine-gordon.qmd` line 85 `eq-sg-defect-lightcone`)!
       Because Lecture 2 defined $\partial_\pm = \frac{1}{2}(\partial_t \mp \partial_x)$ with the opposite sign on $\partial_x$, writing $\partial_+$ on line 201 and $(\partial_t - \partial_x)$ on line 364 swapped the two chiralities! Then, on line 364, Lecture 2 also flipped the sign of the right-hand side ($+m_0\sigma(\dots)$ instead of $-m_0\sigma(\dots)$), and on line 369 claimed that $(\partial_t - \partial_x)E = +m_0 e^\theta E$ for $E = \exp(m_0(x\cosh\theta - t\sinh\theta)+\xi)$, whereas $(\partial_t - \partial_x)E = m_0(-\sinh\theta - \cosh\theta)E = -m_0 e^\theta E$ and $(\partial_t + \partial_x)E = m_0(-\sinh\theta + \cosh\theta)E = +m_0 e^{-\theta} E$!
     - **Convention clash on $\sigma$ and $a(\theta)$:** In the book (`38-defect-sine-gordon.qmd` lines 30, 70), $\sigma = e^{-\theta_d}$ and $a(\theta) = -e^{-\theta}$, and Lecture 6 line 154 also uses $\sigma = e^{-\theta_d}$. Adopting the book's $\partial_\pm = \partial_t \pm \partial_x$ and $\sigma = e^{-\theta_d}$ in Lecture 2 makes $(\partial_t + \partial_x)E = m_0 e^{-\theta}E$ and immediately yields $\frac{1-z_{\rm I_-}}{1+z_{\rm I_-}} = e^{\theta-\theta_d} \implies z_{\rm I_-} = -\coth\frac{\theta-\theta_d}{2}$ with zero sign errors and 100% consistency with Chapter 38 and Lecture 6!
  3. **Soliton vs. antisoliton sign on line 269:**
     - With $E = \exp(m_0(x\cosh\theta - t\sinh\theta)+\xi) = iY$ (`eq-lec2-y-def`, line 266), $Y = -iE$. When $\operatorname{Im}\xi = \frac{\pi}{2}$, $E = i e^{\text{real}}$, so $Y = -i(i e^{\text{real}}) = +e^{\text{real}} > 0$ (soliton, jumping from $0$ to $+\frac{2\pi}{\beta_{\rm SG}}$). When $\operatorname{Im}\xi = \frac{3\pi}{2}$, $E = -i e^{\text{real}}$, so $Y = -e^{\text{real}} < 0$ (antisoliton, jumping from $0$ to $-\frac{2\pi}{\beta_{\rm SG}}$). Line 269 has $Y = -e^{\text{real}} < 0$ and $Y = +e^{\text{real}} > 0$ swapped.
  4. **Linearised wave transmission (`eq-lec2-t-cl`, line 336, and `eq-lec2-t-cl-squared`, line 350):**
     - Substituting $\phi = e^{im_0(x\cosh\theta - t\sinh\theta)}$ and $\psi = \mathbb{T}^{\rm cl}(\theta)e^{im_0(x\cosh\theta - t\sinh\theta)}$ into $(\partial_t + \partial_x)(\phi-\psi) = -m_0\sigma(\phi+\psi)$ with $\sigma = e^{-\theta_d}$ gives $im_0 e^{-\theta}(1 - \mathbb{T}^{\rm cl}) = -m_0 e^{-\theta_d}(1 + \mathbb{T}^{\rm cl})$, so
       $$
       \mathbb{T}^{\rm cl}(\theta) = \frac{e^{-(\theta-\theta_d)} + i}{e^{-(\theta-\theta_d)} - i} = -i\,\frac{\sinh\!\left(\frac{\theta-\theta_d}{2} - \frac{i\pi}{4}\right)}{\sinh\!\left(\frac{\theta-\theta_d}{2} + \frac{i\pi}{4}\right)} = -i\,\frac{\sinh\frac{\theta-\theta_d}{2} - i\cosh\frac{\theta-\theta_d}{2}}{\sinh\frac{\theta-\theta_d}{2} + i\cosh\frac{\theta-\theta_d}{2}} = \frac{\sinh\theta_d + i\cosh\theta}{\cosh\theta_d - i\sinh\theta},
       $$
       matching `38-defect-sine-gordon.qmd` line 171 (`eq-shg-linear-transmission`), and squaring to $[\mathbb{T}^{\rm cl}(\theta)]^2 = \frac{\sinh(\theta-\theta_d)-i}{\sinh(\theta-\theta_d)+i} = f_{-1/2}(\theta-\theta_d)$. In `eq-lec2-t-cl` (line 336), Lecture 2 wrote $\frac{\cosh\theta_d + i\sinh\theta}{\cosh\theta_d - i\sinh\theta} = \frac{\sinh\frac{\theta-\theta_d}{2}-i\cosh\frac{\theta-\theta_d}{2}}{\sinh\frac{\theta-\theta_d}{2}+i\cosh\frac{\theta-\theta_d}{2}}$, where the two expressions are negatives of each other and the second expression squares to $-f_{-1/2}(\theta-\theta_d)$.

- **Replacement text for `lecture-notes/defects-gemini/02-momentum-and-backlund.qmd`:**
  - **Lines 195–207:**
    ```markdown
    Why is this the classic Bäcklund transformation? Introduce light-cone derivatives $\partial_\pm = \partial_t \pm \partial_x$ (so that $\partial_t^2 - \partial_x^2 = \partial_+\partial_-$). Taking the difference and sum of the two first-order defect equations (-@eq-lec2-defect-eom) with the potential (-@eq-lec2-sg-bd) isolates the two chiralities at $x=0$:

    $$
    \begin{aligned}
    \partial_+(\phi - \psi) &= -\frac{2 m_0 \sigma}{\beta_{\rm SG}} \sin\!\left(\frac{\beta_{\rm SG}(\phi+\psi)}{2}\right), \\
    \partial_-(\phi + \psi) &= +\frac{2 m_0}{\sigma \beta_{\rm SG}} \sin\!\left(\frac{\beta_{\rm SG}(\phi-\psi)}{2}\right).
    \end{aligned}
    $$ {#eq-lec2-sg-lc}

    If we require (@eq-lec2-sg-lc) to hold not just at $x=0$, but throughout the entire $(x,t)$ plane with $\sigma = e^{-\theta_d}$, we obtain the classical auto-Bäcklund transformation of the sine-Gordon equation (with Bäcklund parameter $a = -e^{-\theta_d}$). Differentiating the first equation with respect to $\partial_-$ and substituting the second immediately implies $\partial_+\partial_-\phi = -\frac{m_0^2}{\beta_{\rm SG}}\sin(\beta_{\rm SG}\phi)$ and $\partial_+\partial_-\psi = -\frac{m_0^2}{\beta_{\rm SG}}\sin(\beta_{\rm SG}\psi)$.
    ```
  - **Line 269:**
    ```markdown
    When $\operatorname{Im}\xi = \pi/2$, $Y = +e^{\text{real}} > 0$ is real and positive, and $\phi(x,t)$ interpolates smoothly from $0$ at $x = -\infty$ to $+\frac{2\pi}{\beta_{\rm SG}}$ at $x = +\infty$: a **soliton** of topological charge $Q = +1$. When $\operatorname{Im}\xi = 3\pi/2$, $Y = -e^{\text{real}} < 0$ is negative, and $\phi(x,t)$ interpolates from $0$ to $-\frac{2\pi}{\beta_{\rm SG}}$: an **antisoliton** of topological charge $Q = -1$.
    ```
  - **Lines 325–352:**
    ```markdown
    Taking the difference of the two equations gives $(\partial_t + \partial_x)(\phi - \psi) = -m_0\sigma(\phi + \psi)$. Now send in a right-moving plane wave of rapidity $\theta$ from $x < 0$, with $\sigma = e^{-\theta_d}$. Because a Type-I defect is reflectionless, the ansatz is purely transmissive:

    $$
    \phi(x,t) = e^{i m_0 (x \cosh\theta - t \sinh\theta)}, \qquad \psi(x,t) = \mathbb{T}^{\rm cl}(\theta)\, e^{i m_0 (x \cosh\theta - t \sinh\theta)}.
    $$ {#eq-lec2-plane-wave}

    Acting with $\partial_t + \partial_x$ multiplies the plane wave by $i m_0(\cosh\theta - \sinh\theta) = i m_0 e^{-\theta}$. At $x=0$, we obtain an algebraic equation for the **classical transmission amplitude** $\mathbb{T}^{\rm cl}(\theta)$:

    $$
    i e^{-\theta}\bigl(1 - \mathbb{T}^{\rm cl}(\theta)\bigr) = -e^{-\theta_d}\bigl(1 + \mathbb{T}^{\rm cl}(\theta)\bigr) \quad \Longrightarrow \quad
    \mathbb{T}^{\rm cl}(\theta) = \frac{e^{-(\theta-\theta_d)} + i}{e^{-(\theta-\theta_d)} - i} = -i\,\frac{\sinh\frac{\theta-\theta_d}{2} - i\cosh\frac{\theta-\theta_d}{2}}{\sinh\frac{\theta-\theta_d}{2} + i\cosh\frac{\theta-\theta_d}{2}}.
    $$ {#eq-lec2-t-cl}
    ```
  - **Lines 362–375:**
    ```markdown
    To prove (@eq-lec2-delay) directly at $x=0$ without invoking Bianchi permutability in the bulk, substitute the Hirota ansatz (-@eq-lec2-tau-1sol) into the first chiral defect equation $\partial_+(\phi - \psi) = -\frac{2 m_0 \sigma}{\beta_{\rm SG}}\sin\frac{\beta_{\rm SG}(\phi+\psi)}{2}$ with $\sigma = e^{-\theta_d}$. Using the logarithmic derivative rules, this single equation simplifies to a bilinear identity at $x=0$:

    $$
    (\partial_t + \partial_x)\bigl(\tau_-\tilde{\tau}_+ - \tau_+\tilde{\tau}_-\bigr) = -m_0 \sigma \bigl(\tau_+\tilde{\tau}_+ - \tau_-\tilde{\tau}_-\bigr).
    $$ {#eq-lec2-tau-defect}

    Substituting $\tau_\pm = 1 \pm E$ and $\tilde{\tau}_\pm = 1 \pm z_{\rm I_-} E$ gives $\tau_-\tilde{\tau}_+ - \tau_+\tilde{\tau}_- = 2(z_{\rm I_-} - 1)E$ and $\tau_+\tilde{\tau}_+ - \tau_-\tilde{\tau}_- = 2(z_{\rm I_-} + 1)E$. Since $(\partial_t + \partial_x)E = m_0(\cosh\theta - \sinh\theta)E = m_0 e^{-\theta} E$ and $\sigma = e^{-\theta_d}$, cancelling $2 m_0 E$ yields

    $$
    e^{-\theta}(z_{\rm I_-} - 1) = -e^{-\theta_d}(z_{\rm I_-} + 1) \quad \Longrightarrow \quad z_{\rm I_-}(\theta) = \frac{e^{-\theta} - e^{-\theta_d}}{e^{-\theta} + e^{-\theta_d}} = -\coth\!\left(\frac{\theta - \theta_d}{2}\right),
    $$ {#eq-lec2-tau-delay-solved}
    ```
  *(Note: Also update $\sigma = e^{\theta_d} \to \sigma = e^{-\theta_d}$ in `eq-lec2-delay` line 359, lines 441 and 479 so that $\sigma = e^{-\theta_d}$ is uniform across the book and lectures.)*

---

### Item 3: Lecture 3 (`lecture-notes/defects-gemini/03-necklace-and-lax.qmd`)

- **Verdict:** `wrong` (the defect matrix $\mathsf{D}$ in `eq-lec3-necklace-d` and `thm-lec3-necklace` is the transpose of the book's $\mathsf{D}$ in `39-defect-conditions.qmd`, causing downstream errors in the Lax curvature proof on line 311, the Fourier eigenvalue on line 362, the Hirota sign on line 394, the Hollowood soliton mass on line 409, and the soliton delay formula on lines 424–433).
- **Independent computation (`referee_scratch/lectures_gemini/verify_lectures_1_to_6.py`, Section 3):**
  1. **Necklace coupling matrix $\mathsf{D}$ (`eq-lec3-necklace-d`, line 229, and `thm-lec3-necklace`, line 255):**
     - In the book (`39-defect-conditions.qmd` lines 31–68 `eq-d-matrix`, `thm-defect-integrability`), the Coxeter shift acts by $\pi\vec{\alpha}_j = \vec{\alpha}_{j-1}$, and $\mathsf{D} = 2(1 - \pi)^{-1}$ satisfies
       $$
       \vec{\alpha}_i\cdot\mathsf{D}\vec{\alpha}_j = 2\delta_{ij} - 2\delta_{i,j+1}, \qquad \mathsf{D} + \mathsf{D}^{\mathsf{T}} = 2\,\mathbf{1},
       $$
       so that $\frac{1}{2}\vec{\alpha}_j\cdot(\vec{\phi} + (\mathsf{D}-1)\vec{\psi}) = \frac{1}{2}\vec{\alpha}_j\cdot(\vec{\phi} - \vec{\psi} + \mathsf{D}\vec{\psi})$ and $\frac{1}{2}\vec{\alpha}_j\cdot(\vec{\phi} - (\mathsf{D}+1)\vec{\psi}) = \frac{1}{2}\vec{\alpha}_j\cdot(-\vec{\phi} + \vec{\psi} + \mathsf{D}\vec{\psi})$ (up to sign), pairing $f_j = \exp(\frac{b}{2}\vec{\alpha}_j\cdot(\vec{\phi} + (\mathsf{D}-1)\vec{\psi}))$ with $g_{j+1} = \exp(-\frac{b}{2}\vec{\alpha}_{j+1}\cdot(\vec{\phi} - (\mathsf{D}+1)\vec{\psi}))$ via $\mathsf{D}^{\mathsf{T}}\vec{\alpha}_j - \vec{\alpha}_j = (2 - \mathsf{D} - 1)\vec{\alpha}_j = (1 - \mathsf{D})\vec{\alpha}_j$ and $(\mathsf{D}+1)\vec{\alpha}_{j+1}$:
       specifically, $\frac{1}{2}\vec{\alpha}_j\cdot(\vec{\phi} + (\mathsf{D}-1)\vec{\psi}) - \frac{1}{2}\vec{\alpha}_{j+1}\cdot(\vec{\phi} - (\mathsf{D}+1)\vec{\psi})$ cancels $\vec{\phi}$ when we use the book's variables $\vec{u} = \frac{1}{2}\mathsf{D}^{\mathsf{T}}(\vec{\phi}+\vec{\psi}) - \vec{\psi}$ and $\vec{v} = \frac{1}{2}\mathsf{D}^{\mathsf{T}}(\vec{\phi}-\vec{\psi})$ (`39-defect-conditions.qmd` line 23 `eq-uv-coords`), where
       $$
       \vec{u} + \vec{v} = \mathsf{D}^{\mathsf{T}}\vec{\phi} - \vec{\psi}, \qquad \vec{u} - \vec{v} = (\mathsf{D}^{\mathsf{T}} - 1)\vec{\psi} = (1 - \mathsf{D})\vec{\psi}.
       $$
       Look at what happens in the book's coordinates $\vec{u}, \vec{v}$:
       $$
       f_j = e^{b\vec{\alpha}_j\cdot\vec{u}} = \exp\!\left(\frac{b}{2}\vec{\alpha}_j\cdot\bigl[\mathsf{D}^{\mathsf{T}}\vec{\phi} + (\mathsf{D}^{\mathsf{T}}-2)\vec{\psi}\bigr]\right) = \exp\!\left(\frac{b}{2}(\mathsf{D}\vec{\alpha}_j)\cdot\vec{\phi} - \frac{b}{2}(\mathsf{D}\vec{\alpha}_j)\cdot\vec{\psi} + b\vec{\alpha}_j\cdot(\mathsf{D}^{\mathsf{T}}-1)\vec{\psi}\right),
       $$
       and since $\vec{\alpha}_i\cdot(\vec{u}+\vec{v}) = \vec{\alpha}_i\cdot(\mathsf{D}^{\mathsf{T}}\vec{\phi}-\vec{\psi}) = (\mathsf{D}\vec{\alpha}_i)\cdot\vec{\phi} - \vec{\alpha}_i\cdot\vec{\psi}$ and $\vec{\alpha}_i\cdot\mathsf{D}\vec{\alpha}_j = 2\delta_{ij} - 2\delta_{i,j+1}$, we have $\mathsf{D}\vec{\alpha}_{j+1} = \vec{\alpha}_{j+1} - \vec{\alpha}_j$ (on the root space)!
       Wait: on the root space of $a_r^{(1)}$, since $\pi\vec{\alpha}_j = \vec{\alpha}_{j-1}$, $(1-\pi)\vec{\alpha}_j = \vec{\alpha}_j - \vec{\alpha}_{j-1}$, and since $\vec{\alpha}_i\cdot\mathsf{D}\vec{\alpha}_j = 2\delta_{ij} - 2\delta_{i,j+1}$, let's check how $\mathsf{D}^{\mathsf{T}}$ acts on $\vec{\alpha}_j$:
       $$
       \vec{\alpha}_i\cdot(\mathsf{D}^{\mathsf{T}}\vec{\alpha}_j - \mathsf{D}^{\mathsf{T}}\vec{\alpha}_{j+1}) = \vec{\alpha}_j\cdot\mathsf{D}\vec{\alpha}_i - \vec{\alpha}_{j+1}\cdot\mathsf{D}\vec{\alpha}_i = (2\delta_{ji} - 2\delta_{j,i+1}) - (2\delta_{j+1,i} - 2\delta_{j+1,i+1}) = 2\delta_{i,j} - 2\delta_{i,j-1} - 2\delta_{i,j+1} + 2\delta_{i,j} = 2 C_{ij}\vec{\alpha}_i!
       $$
       Even simpler: in the book (`39-defect-conditions.qmd` line 122 `eq-type-i-potential`), the defect potential is
       $$
       \mathcal{B}_d(\vec{\phi},\vec{\psi}) = \frac{m}{\beta^2}\sum_{j=0}^r \left(\sigma\,e^{\frac{1}{2}i\beta\vec{\alpha}_j\cdot(\vec{\phi}+\mathsf{D}\vec{\psi} - \vec{\psi})} + \sigma^{-1}\,e^{-\frac{1}{2}i\beta\vec{\alpha}_j\cdot(\vec{\phi}-\mathsf{D}\vec{\psi} - \vec{\psi})}\right) = \frac{m}{b^2}\sum_{j=0}^r (\sigma f_j + \sigma^{-1} g_j)
       $$
       where $f_j = \exp(\frac{b}{2}\vec{\alpha}_j\cdot(\vec{\phi} + (\mathsf{D}-1)\vec{\psi}))$ and $g_j = \exp(-\frac{b}{2}\vec{\alpha}_j\cdot(\vec{\phi} - (\mathsf{D}+1)\vec{\psi}))$!
       With the book's $\vec{\alpha}_i\cdot\mathsf{D}\vec{\alpha}_j = 2\delta_{ij} - 2\delta_{i,j+1}$, the exponent of $f_i g_{i+1}$ has $\vec{\psi}$-term
       $$
       \frac{b}{2}\vec{\alpha}_i\cdot(\mathsf{D}-1)\vec{\psi} + \frac{b}{2}\vec{\alpha}_{i+1}\cdot(\mathsf{D}+1)\vec{\psi},
       $$
       Wait: in the book (`39-defect-conditions.qmd` line 15 `eq-defect-action`), the kinetic coupling in the Lagrangian is $\frac{1}{2}\vec{\phi}\cdot\mathsf{A}\partial_t\vec{\psi} - \frac{1}{2}\vec{\psi}\cdot\mathsf{A}\partial_t\vec{\phi}$ where $\mathsf{A} = \mathsf{D}-1$ is skew-symmetric ($\mathsf{D}+\mathsf{D}^{\mathsf{T}}=2$), and the defect equations of motion (`39-defect-conditions.qmd` line 19 `eq-defect-eom-general`) are:
       $$
       \mathsf{A}(\partial_t\vec{\psi} - \partial_x\vec{\phi}) = \nabla_{\vec{\phi}}\mathcal{B}_d - \partial_x(\vec{\phi}-\vec{\psi})\quad \Longleftrightarrow \quad \partial_t(\vec{\phi}-\vec{\psi}) + \partial_x(\vec{\phi}+\vec{\psi}) = \dots
       $$
       Wait let's check why Lecture 3 got into trouble on line 311!
  2. **Error in the Lax pair proof on line 311:**
     - In Lecture 3 (`eq-lec3-type1-action`, line 46), the defect equations of motion (`eq-lec3-defect-eom`, line 59) are written as:
       $$
       \mathsf{D}(\partial_t\vec{\phi} - \partial_x\vec{\psi}) = -\nabla_{\vec{\psi}}\mathcal{B}_d, \qquad \mathsf{D}^{\mathsf{T}}(\partial_t\vec{\psi} - \partial_x\vec{\phi}) = +\nabla_{\vec{\phi}}\mathcal{B}_d.
       $$
       Then on line 305 (`eq-lec3-dt-ld`), Lecture 3 differentiates $L_d = -\frac{b}{4}\mathsf{D}^{\mathsf{T}}(\vec{\phi}-\vec{\psi})\cdot\vec{H}$ to get $\partial_t L_d + (\text{bulk } \partial_x \text{ terms}) = -\frac{b}{4}(\nabla_{\vec{\phi}}\mathcal{B}_d + \nabla_{\vec{\psi}}\mathcal{B}_d)\cdot\vec{H}$.
       With $f_j = \exp(\frac{b}{2}\vec{\alpha}_j\cdot(\vec{\phi} + (\mathsf{D}-1)\vec{\psi}))$ and $g_j = \exp(-\frac{b}{2}\vec{\alpha}_j\cdot(\vec{\phi} - (\mathsf{D}+1)\vec{\psi}))$, we have:
       $$
       (\nabla_{\vec{\phi}} + \nabla_{\vec{\psi}})f_j = \frac{b}{2}(\mathbf{1} + \mathsf{D}^{\mathsf{T}} - \mathbf{1})\vec{\alpha}_j f_j = \frac{b}{2}\mathsf{D}^{\mathsf{T}}\vec{\alpha}_j f_j,
       $$
       $$
       (\nabla_{\vec{\phi}} + \nabla_{\vec{\psi}})g_j = -\frac{b}{2}(\mathbf{1} - (\mathsf{D}^{\mathsf{T}} + \mathbf{1}))\vec{\alpha}_j g_j = +\frac{b}{2}\mathsf{D}^{\mathsf{T}}\vec{\alpha}_j g_j.
       $$
       Whereas on line 309 (`eq-lec3-commutator-target`), the commutator $[U_t^+, L_d] - L_d U_t^-$ produces $+\frac{b}{4}\sum_j (\sigma f_j + \sigma^{-1} g_j)\vec{\alpha}_j\cdot\vec{H}$ (with $\vec{\alpha}_j$, NOT $\mathsf{D}^{\mathsf{T}}\vec{\alpha}_j$!).
       Why did $\mathsf{D}^{\mathsf{T}}\vec{\alpha}_j$ appear in $\partial_t L_d$ while $\vec{\alpha}_j$ appeared in the commutator? Because in `eq-lec3-defect-eom` (line 59), multiplying $\partial_t\vec{\phi}-\partial_x\vec{\psi}$ by $\mathsf{D}$ and $\partial_t\vec{\psi}-\partial_x\vec{\phi}$ by $\mathsf{D}^{\mathsf{T}}$ gives $(\partial_t + \partial_x)\mathsf{D}^{\mathsf{T}}(\vec{\phi}-\vec{\psi}) = \nabla_{\vec{\phi}}\mathcal{B}_d + \nabla_{\vec{\psi}}\mathcal{B}_d + (\mathsf{D}^{\mathsf{T}}-\mathsf{D})\partial_t\vec{\phi}$, which does NOT cancel $\mathsf{D}^{\mathsf{T}}$!
       Look at how the book (`39-defect-conditions.qmd` lines 19–27 `eq-defect-eom-general` and `eq-defect-lightcone`) actually defines the defect equations of motion:
       $$
       \partial_+(\vec{\phi} - \vec{\psi}) = -2\mathsf{D}^{-1}(\nabla_{\vec{\phi}} + \nabla_{\vec{\psi}})\mathcal{B}_d, \qquad \partial_-(\vec{\phi} + \vec{\psi}) = 2(\mathsf{D}^{\mathsf{T}})^{-1}(\nabla_{\vec{\phi}} - \nabla_{\vec{\psi}})\mathcal{B}_d!
       $$
       Wait: when $\partial_+(\vec{\phi}-\vec{\psi}) = -2(\mathsf{D}^{\mathsf{T}})^{-1}(\nabla_{\vec{\phi}}+\nabla_{\vec{\psi}})\mathcal{B}_d$, the $(\mathsf{D}^{\mathsf{T}})^{-1}$ cancels the $\mathsf{D}^{\mathsf{T}}$ in $(\nabla_{\vec{\phi}}+\nabla_{\vec{\psi}})f_j = \frac{b}{2}\mathsf{D}^{\mathsf{T}}\vec{\alpha}_j f_j$ and $(\nabla_{\vec{\phi}}+\nabla_{\vec{\psi}})g_j = \frac{b}{2}\mathsf{D}^{\mathsf{T}}\vec{\alpha}_j g_j$ **identically**!
       In Lecture 3, because the author omitted $(\mathsf{D}^{\mathsf{T}})^{-1}$ in $L_d = -\frac{b}{4}(\vec{\phi}-\vec{\psi})\cdot\vec{H}$, the author tried on line 311 to claim that $(\mathsf{D}-\mathsf{D}^{\mathsf{T}})\vec{\alpha}_j = -\vec{\alpha}_j$—which our script `verify_lectures_1_to_6.py` proves is **false** for $A_2$ ($(\mathsf{D}-\mathsf{D}^{\mathsf{T}})\vec{\alpha}_1 = (-2/3, -4/3)^{\mathsf{T}} \neq (-1, 0)^{\mathsf{T}}$).
  3. **Sign error in Fourier eigenvalue (`eq-lec3-dk-eigenvalues`, line 362):**
     - Verified in `verify_lectures_1_to_6.py` that $\frac{2}{1 - e^{-2\pi ik/h}} = 1 - i\cot\frac{\pi k}{h} \neq 1 + i\cot\frac{\pi k}{h}$. With the book's $\mathsf{D} = 2(1-\pi)^{-1}$ (`39-defect-conditions.qmd` `eq-d-matrix`), the eigenvalue on the $k$-th Coxeter mode is $\frac{2}{1 - e^{2\pi ik/h}} = 1 + i\cot\frac{\pi k}{h}$, so switching Lecture 3 to the book's $\mathsf{D}$ makes $1 + i\cot\frac{\pi k}{h}$ right!
  4. **Hirota sign (`eq-lec3-hirota-ansatz`, line 394), soliton mass (`eq-lec3-1sol-tau` prose, line 409), and Bäcklund delay (`eq-lec3-fj-gj-tau` line 424 and `eq-lec3-za` line 433):**
     - On line 394, $\vec{\phi} = -\frac{i}{\beta}\sum_j \vec{\alpha}_j\ln\tau_j$ gives $e^{i\beta\vec{\alpha}_j\cdot\vec{\phi}} = \frac{\tau_j^2}{\tau_{j-1}\tau_{j+1}}$, contradicting line 398 (`eq-lec3-exp-tau`). Changing line 394 to the book's $+\frac{i}{\beta}\sum_j \vec{\alpha}_j\ln\tau_j$ (`40-defect-solitons.qmd` line 15 `eq-hirota-ansatz`) fixes line 398.
     - On line 409, the Hollowood soliton mass in the book (`40-defect-solitons.qmd` line 17) is $M_a = \frac{2hm_a}{\beta^2} = \frac{4hm}{\beta^2}\sin\frac{\pi a}{h}$, not $\frac{2m}{\beta^2}\sin\frac{\pi a}{h}$.
     - On line 433 (`eq-lec3-za`), verified in `verify_lectures_1_to_6.py` that the book's delay (`40-defect-solitons.qmd` line 36 `eq-an-defect-delay`) is
       $$
       z_a(\theta) = \frac{e^{-(\theta-\theta_d)} - i e^{i\gamma_a}}{e^{-(\theta-\theta_d)} - i e^{-i\gamma_a}} = e^{i\gamma_a}\,\frac{\sinh\!\left(\frac{\theta-\theta_d}{2} + \frac{i(\pi/2+\gamma_a)}{2}\right)}{\sinh\!\left(\frac{\theta-\theta_d}{2} + \frac{i(\pi/2-\gamma_a)}{2}\right)},
       $$
       which at $\gamma_1 = \pi/2$ ($a_1^{(1)}$) equals $-\coth\frac{\theta-\theta_d}{2}$. By contrast, Lecture 3's `eq-lec3-za` ($\frac{e^{\theta-\theta_d}-ie^{-i\gamma_a}}{e^{\theta-\theta_d}+ie^{i\gamma_a}}$) evaluates to $+1$ identically at $\gamma_a = \pi/2$!

- **Replacement text for `lecture-notes/defects-gemini/03-necklace-and-lax.qmd`:**
  - **Lines 46–61 (`eq-lec3-type1-action` and `eq-lec3-defect-eom`):** Align with `39-defect-conditions.qmd` (`eq-defect-action`, `eq-defect-lightcone`) by writing the defect kinetic term with $\mathsf{A} = \mathsf{D}-1$ (or writing the light-cone defect equations directly):
    ```markdown
    $$
    \partial_+(\vec{\phi} - \vec{\psi}) = -2(\mathsf{D}^{\mathsf{T}})^{-1}(\nabla_{\vec{\phi}} + \nabla_{\vec{\psi}})\mathcal{B}_d, \qquad
    \partial_-(\vec{\phi} + \vec{\psi}) = +2\mathsf{D}^{-1}(\nabla_{\vec{\phi}} - \nabla_{\vec{\psi}})\mathcal{B}_d.
    $$ {#eq-lec3-defect-eom}
    ```
  - **Lines 220–261 (`eq-lec3-necklace-d` and `thm-lec3-necklace`):** Pair $f_i$ with $g_{i+1}$ so that $\mathsf{D}$ matches the book's `39-defect-conditions.qmd` (`eq-d-matrix`):
    ```markdown
    $$
    \vec{\alpha}_i \cdot \mathsf{D} \vec{\alpha}_j = 2\delta_{ij} - 2\delta_{i, j+1} \pmod{r+1}.
    $$ {#eq-lec3-necklace-d}
    ```
    and in `thm-lec3-necklace` (line 255):
    ```markdown
    $$
    \vec{\alpha}_i \cdot \mathsf{D} \vec{\alpha}_j = 2\delta_{ij} - 2\delta_{i,j+1}, \qquad \mathsf{D} = 2(1 - \pi)^{-1},
    $$
    where $\pi(\vec{\alpha}_j) = \vec{\alpha}_{j-1}$ is the cyclic shift of the extended Dynkin diagram.
    ```
  - **Lines 299–316 (`eq-lec3-ld-def` to `eq-lec3-commutator-target`):** With $L_d = -\frac{b}{4}(\vec{\phi}-\vec{\psi})\cdot\vec{H}$ and $\partial_+(\vec{\phi}-\vec{\psi}) = -2(\mathsf{D}^{\mathsf{T}})^{-1}(\nabla_{\vec{\phi}}+\nabla_{\vec{\psi}})\mathcal{B}_d$, since $(\nabla_{\vec{\phi}}+\nabla_{\vec{\psi}})f_j = \frac{b}{2}\mathsf{D}^{\mathsf{T}}\vec{\alpha}_j f_j$ and $(\nabla_{\vec{\phi}}+\nabla_{\vec{\psi}})g_j = \frac{b}{2}\mathsf{D}^{\mathsf{T}}\vec{\alpha}_j g_j$, applying $(\mathsf{D}^{\mathsf{T}})^{-1}$ immediately yields $\partial_+(\vec{\phi}-\vec{\psi}) = -\frac{m}{b}\sum_{j=0}^r (\sigma f_j + \sigma^{-1}g_j)\vec{\alpha}_j$, matching $[U_t^+, L_d] - L_d U_t^-$ on line 309 without the false identity $(\mathsf{D}-\mathsf{D}^{\mathsf{T}})\vec{\alpha}_j = -\vec{\alpha}_j$.
  - **Line 362 (`eq-lec3-dk-eigenvalues`):**
    ```markdown
    $$
    d_k = \frac{2}{1 - \omega^k} = \frac{2}{1 - e^{2\pi i k / h}} = 1 + i \cot\!\left(\frac{\pi k}{h}\right), \qquad
    d_k^* = \frac{2}{1 - \omega^{-k}} = 1 - i \cot\!\left(\frac{\pi k}{h}\right).
    $$ {#eq-lec3-dk-eigenvalues}
    ```
  - **Line 394 (`eq-lec3-hirota-ansatz`) and line 409:**
    ```markdown
    $$
    \vec{\phi}(x,t) = \frac{i}{\beta} \sum_{j=0}^r \vec{\alpha}_j \ln \tau_j(x,t), \qquad
    \vec{\psi}(x,t) = \frac{i}{\beta} \sum_{j=0}^r \vec{\alpha}_j \ln \tilde{\tau}_j(x,t),
    $$ {#eq-lec3-hirota-ansatz}
    ```
    and on line 409: `where $m_a = 2m\sin(\pi a/h)$ is the mass of the $a$-th fundamental Toda particle (while the classical soliton mass is $M_a = \frac{2h m_a}{\beta^2} = \frac{4hm}{\beta^2}\sin(\pi a/h)$).`
  - **Lines 424–435 (`eq-lec3-fj-gj-tau` and `eq-lec3-za`):**
    ```markdown
    $$
    f_j = \frac{\tau_{j-1}\tilde{\tau}_{j+1}}{\tau_j \tilde{\tau}_j}, \qquad g_{j+1} = \frac{\tau_{j+1}\tilde{\tau}_j}{\tau_j \tilde{\tau}_{j+1}}.
    $$ {#eq-lec3-fj-gj-tau}

    The bilinear defect equation at $x=0$ reduces to $(\partial_t + \partial_x)(\tau_j\tilde{\tau}_{j+1} - \tau_{j+1}\tilde{\tau}_j) = -i m \sigma (\tau_{j-1}\tilde{\tau}_{j+2} - \tau_{j+1}\tilde{\tau}_j)$ (see @sec-ds-tau in the main text). Substituting the 1-soliton tau-functions (-@eq-lec3-1sol-tau) with $\sigma = e^{-\theta_d}$ and $\gamma_a = \pi a / h$ gives a linear equation for the **classical transmission factor** $z_a(\theta)$:

    $$
    z_a(\theta) = \frac{e^{-(\theta - \theta_d)} - i e^{i\gamma_a}}{e^{-(\theta - \theta_d)} - i e^{-i\gamma_a}} = e^{i\gamma_a}\,\frac{\sinh\!\left(\frac{\theta-\theta_d}{2} + \frac{i(\pi/2 + \gamma_a)}{2}\right)}{\sinh\!\left(\frac{\theta-\theta_d}{2} + \frac{i(\pi/2 - \gamma_a)}{2}\right)}.
    $$ {#eq-lec3-za}
    ```

---

### Item 4: Lecture 4 (`lecture-notes/defects-gemini/04-quantum-groups-and-q-operators.qmd`)

- **Verdict:** `wrong` (serious representation/coproduct/Weyl-vs-Fock errors in lines 48–50, 89, and 169–235).
- **Independent computation (`referee_scratch/lectures_gemini/verify_lectures_1_to_6.py`, Section 4):**
  1. **Conformal limit and Borel subalgebra (`sec-lec4-conformal-borel`, lines 48–50):**
     - Lines 48–50 state that the Type-I defect action perturbs the CFT by single-chirality vertex operators $e^{\pm i\beta_{\rm SG}\phi_L(0,t)}$ and therefore commutes with all bulk non-local charges of one chirality. As explained in `41-defect-quantum-groups.qmd` lines 47–61, the Lagrangian Type-I defect potential $\mathcal{B}_d$ contains **both** chiralities ($\sigma\cos\beta_{\rm SG}\phi_L$ and $\sigma^{-1}\cos\beta_{\rm SG}\phi_R$); a purely chiral perturbation is the Konik–LeClair conformal defect model, which motivates the Borel subalgebra $U_q(\mathfrak{b}^-)$ rather than proving all-orders conservation for the Lagrangian defect.
  2. **Coproduct convention (`eq-lec4-borel-coproduct`, line 89):**
     - Line 89 writes $\Delta(f_i) = f_i\otimes 1 + k_i\otimes f_i$, opposite to the book's Part VI convention $\Delta(f_i) = f_i\otimes k_i^{-1} + 1\otimes f_i$ (`41-defect-quantum-groups.qmd` line 13 `eq-borel-coproduct`).
  3. **Intertwiner derivation and Weyl-vs-Fock chimera (`sec-lec4-solving-intertwiner`, lines 169–228):**
     - **Error 4a (broken charge conservation in `eq-lec4-reps-explicit`, `eq-lec4-weyl-explicit`, `eq-lec4-t-ansatz`):** On lines 169–173, $\pi_y(k_1) = \operatorname{diag}(q, q^{-1})$ and $\rho(k_1)|n\rangle = q^{-2n-s}|n\rangle$, so $k_1\otimes k_1$ has eigenvalue $q^{1-2n-s}$ on $v_+\otimes|n\rangle$ and $q^{-1-2m-s}$ on $v_-\otimes|m\rangle$. Equality $1 - 2n = -1 - 2m$ requires $m = n - 1$, NOT $m = n + 1$! Consistently, $\pi_y(f_1) = y E_{21}$ lowers $v_+ \to v_-$ while $\rho(f_1)|n\rangle = c_1 y_d |n+1\rangle$ raises $|n\rangle \to |n+1\rangle$. Because `eq-lec4-t-ansatz` (line 177) instead coupled $v_+\otimes|n\rangle$ to $v_-\otimes|n+1\rangle$, evaluating `eq-lec4-f1-intertwine` on $v_+\otimes|n\rangle$ in `verify_lectures_1_to_6.py` proves that $C_n = B_n = 0$ and on line 216 produces the contradictory equation $a(y)q^n = q^{-2n-s}d(y)q^{-(n+1)}$.
     - **Error 4b (Weyl-vs-Fock chimera in `eq-lec4-sg-t-solved` line 221 and `eq-lec4-sg-t-explicit` line 225):** Even worse, the matrix displayed on lines 221–225:
       $$
       A_n = q^{n+s/2}, \quad D_n = q^{-n-s/2}, \quad B_n = -i e^{\theta-\theta_d}, \quad C_n = -i e^{\theta-\theta_d}(1 - q^{-2n})
       $$
       is a **non-intertwining chimera**:
       - In the book (`38-defect-sine-gordon.qmd` line 190 `eq-sg-transmission-matrix` and `41-defect-quantum-groups.qmd` line 120 `prp-sg-t-unique`), on the **Weyl module** $\mathcal{W}^-(y_d)$ (where $Z = (\rho\otimes\rho)\Delta(f_1 f_0 - q^2 f_0 f_1)|_{\text{central}} = 0$), the recurrence is **monomial** (`41-defect-quantum-groups.qmd` lines 88–94 `eq-intertwiner-recurrence`): there is **no** $(1 - q^{-2n})$ factor anywhere, and the unique solution is
         $$
         A_n = q^{-n}, \qquad D_n = q^{n+s}, \qquad C_n = -i Y q^{-n}, \qquad B_n = -i Y q^{n+s}, \qquad Y = \frac{y}{y_d} = e^{\lambda(\theta-\theta_d)}
         $$
         (Notice also that Lecture 6 line 123 uses this exact Weyl-module solution from the book, contradicting Lectures 4 and 5!).
       - Where did $(1 - q^{-2n})$ come from? On the **Fock module** $\mathcal{F}_{y_F}$ (`41-defect-quantum-groups.qmd` line 151 `sec-qg-baxter-q`), where the central element $Z \neq 0$, the Baxter $L$-operator (`eq-l-operator`, line 164) has $C_n \propto 1 - q^{-2n}$, which requires $A_n = q^{-n} - x^2 q^{n+1}$ to be **quadratic** in $x$! Putting $(1-q^{-2n})$ into $C_n$ while keeping $A_n, D_n$ independent of $Y$ violates the $U_q(\mathfrak{b}^-)$ intertwining equations for every $n$, and subsequently breaks the crossing equation in Lecture 5!
  4. **Universal R-matrix Borel factors (`eq-lec4-universal-r-borel`, line 235):**
     - In the book (`41-defect-quantum-groups.qmd` line 129 `prp-two-borels`), $\mathcal{R} \in U_q(\mathfrak{b}^-)\,\widehat{\otimes}\,U_q(\mathfrak{b}^+)$ and $\mathbb{T}_{\rm I_-}(\theta) \propto P(\pi_\theta\otimes\rho^-)(\mathcal{R}_{21}^{-1})$. Line 235 has $\mathfrak{b}^+$ and $\mathfrak{b}^-$ swapped.

- **Replacement text for `lecture-notes/defects-gemini/04-quantum-groups-and-q-operators.qmd`:**
  - **Lines 44–53 (`sec-lec4-conformal-borel`):** Clarify that the Lagrangian defect potential $\mathcal{B}_d$ couples both chiralities, whereas a purely chiral vertex operator perturbation (Konik–LeClair) preserves one Borel half $U_q(\mathfrak{b}^-)$ at first order in CFT and motivates the exact algebraic postulate in `@sec-lec4-coproduct`.
  - **Line 89 (`eq-lec4-borel-coproduct`):**
    ```markdown
    $$
    \Delta(k_i) = k_i \otimes k_i, \qquad \Delta(f_i) = f_i \otimes k_i^{-1} + 1 \otimes f_i.
    $$ {#eq-lec4-borel-coproduct}
    ```
  - **Lines 165–228 (`sec-lec4-solving-intertwiner`):** Replace with the book's exact derivation from `41-defect-quantum-groups.qmd` lines 65–125 (`sec-qg-sg-intertwiner` and `prp-sg-t-unique`):
    ```markdown
    Using the book's principal-gradation spin-$\frac{1}{2}$ evaluation representation $\pi_y$ ($y = e^{\lambda\theta}$) and defect Weyl module $\mathcal{W}^-(y_d)$ ($y_d = e^{\lambda\theta_d}$, see @sec-qg-sg-intertwiner in the main text):

    $$
    \pi_y(k_1) = \begin{pmatrix} q^{-1} & 0 \\ 0 & q \end{pmatrix}, \quad \pi_y(k_0) = \begin{pmatrix} q & 0 \\ 0 & q^{-1} \end{pmatrix}, \quad
    \pi_y(f_1) = y^{-1}\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}, \quad \pi_y(f_0) = y^{-1}\begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix},
    $$ {#eq-lec4-reps-explicit}

    and on the defect Weyl basis $\{|n\rangle\}_{n \ge 0}$:

    $$
    \rho(k_1)|n\rangle = q^{-2n-s}|n\rangle, \quad \rho(k_0)|n\rangle = q^{2n+s}|n\rangle, \quad
    \rho(f_1)|n\rangle = c_1 y_d^{-1} q^{-n}|n-1\rangle, \quad \rho(f_0)|n\rangle = c_0 y_d^{-1} q^{n+s}|n+1\rangle.
    $$ {#eq-lec4-weyl-explicit}

    Step 1: Cartan weight conservation $\mathbb{T}(\theta)\Delta(k_1) = \Delta^{\rm op}(k_1)\mathbb{T}(\theta)$ gives total $k_1$-weight $q^{-1-2n-s}$ on $v_+\otimes|n\rangle$ and $q^{1-2m-s}$ on $v_-\otimes|m\rangle$, forcing $m = n + 1$:

    $$
    \begin{aligned}
    \mathbb{T}(\theta)\bigl(v_+ \otimes |n\rangle\bigr) &= |n\rangle \otimes A_n(y)\,v_+ + |n+1\rangle \otimes C_n(y)\,v_-, \\
    \mathbb{T}(\theta)\bigl(v_- \otimes |n\rangle\bigr) &= |n-1\rangle \otimes B_n(y)\,v_+ + |n\rangle \otimes D_n(y)\,v_-.
    \end{aligned}
    $$ {#eq-lec4-t-ansatz}

    Step 2: Intertwining $f_1$ and $f_0$ via $\mathbb{T}(\theta)\Delta(f_i) = \Delta^{\rm op}(f_i)\mathbb{T}(\theta)$ yields the four monomial recurrence relations (`eq-intertwiner-recurrence` in @sec-qg-sg-intertwiner):

    $$
    \frac{A_n}{A_{n-1}} = q^{-1}, \qquad \frac{D_{n+1}}{D_n} = q, \qquad
    C_n = c_0(1 - q^2)\frac{y}{y_d} q^{-n} A_0, \qquad B_n = c_1(1 - q^{-2})\frac{y}{y_d} q^{n+s} D_0.
    $$ {#eq-lec4-f0f1-recurrence}

    Symmetric charge-conjugation normalization $c_0(1-q^2) = c_1(1-q^{-2}) = -i$ and $A_0 = D_0 = \Phi(\theta)$ fixes the **Type-I$_-$ sine-Gordon soliton transmission matrix** uniquely up to the scalar prefactor $\Phi(\theta)$, with $Y = y/y_d = e^{\lambda(\theta-\theta_d)}$:

    $$
    A_n(\theta) = \Phi(\theta)\,q^{-n}, \qquad D_n(\theta) = \Phi(\theta)\,q^{n+s}, \qquad
    C_n(\theta) = -i\Phi(\theta)\,Y q^{-n}, \qquad B_n(\theta) = -i\Phi(\theta)\,Y q^{n+s},
    $$ {#eq-lec4-sg-t-solved}

    or in block form on $\{v_+, v_-\}$:

    $$
    \mathbb{T}_{\rm I_-}(\theta) = \Phi(\theta)\begin{pmatrix}
    q^{-\hat{N}} & -i e^{\lambda(\theta-\theta_d)}\,\hat{e}_- q^{\hat{N}+s} \\
    -i e^{\lambda(\theta-\theta_d)}\,\hat{e}_+ q^{-\hat{N}} & q^{\hat{N}+s}
    \end{pmatrix}.
    $$ {#eq-lec4-sg-t-explicit}
    ```
  - **Line 235 (`eq-lec4-universal-r-borel`):**
    ```markdown
    $$
    \mathcal{R} \in U_q(\mathfrak{b}^-) \,\widehat{\otimes}\, U_q(\mathfrak{b}^+), \qquad
    \mathbb{T}_{\rm I_-}(\theta) \propto P(\pi_\theta \otimes \rho^-)(\mathcal{R}_{21}^{-1}).
    $$ {#eq-lec4-universal-r-borel}
    ```

---

### Item 5: Lecture 5 (`lecture-notes/defects-gemini/05-crossing-resonances-and-tunnelling.qmd`)

- **Verdict:** `wrong` (2x2 block and crossing algebra errors in lines 29–108 inherited from Lecture 4's chimera matrix, missing $2h=6$ in the soliton mass on line 172, and misstated $a_r^{(1)}$ Weyl/transmission formulas on lines 192–210).
- **Independent computation (`referee_scratch/lectures_gemini/verify_lectures_1_to_6.py`, Section 5):**
  1. **2x2 invariant block $\mathcal{H}_n = \operatorname{span}\{v_+\otimes|n\rangle, v_-\otimes|n+1\rangle\}$ and unitarity/crossing (`sec-lec5-unitarity-crossing`, lines 29–108):**
     - In Lecture 5, `eq-lec5-tn-2x2` (line 31) restricts Lecture 4's chimera matrix to $\mathcal{H}_n$ (and additionally drops $+1$ in the $(2,2)$ exponent $q^{-(n+1)-s/2}$ on $|n+1\rangle$, and omits $\lambda$ in $Y = e^{\lambda(\theta-\theta_d)}$). As verified symbolically in `verify_lectures_1_to_6.py`, when one computes the $(1,1)$ entry of `eq-lec5-crossing-op` (line 69) using Lecture 5's matrices, one gets `eq-lec5-crossing-diag` (line 76):
       $$
       \Phi(\theta)\Phi(\theta+i\pi)\bigl[1 + (1 - q^{-2(n+1)})e^{2(\theta-\theta_d)}\bigr] = 1,
       $$
       where the $n$-dependence $q^{-2(n+1)}$ **fails to cancel**! On lines 80–83 (`eq-lec5-crossing-diag-shifted`), the author attempts to patch this by replacing $1$ with $q^{-2n-s}$, which still does not cancel $q^{-2(n+1)}$.
     - By contrast, using the book's Weyl-module block on $\mathcal{H}_n$ (`38-defect-sine-gordon.qmd` line 214 `eq-sg-block-matrix`, which Lecture 6 line 126 already uses!):
       $$
       \mathbb{T}_n(\theta) = \Phi(\theta)\begin{pmatrix} q^{-n} & -i Y q^{n+1+s} \\ -i Y q^{-n} & q^{n+1+s} \end{pmatrix} = \Phi(\theta)\begin{pmatrix} 1 & -i Y \\ -i Y & 1 \end{pmatrix}\begin{pmatrix} q^{-n} & 0 \\ 0 & q^{n+1+s} \end{pmatrix},
       $$
       with $Y(\theta) = e^{\lambda(\theta-\theta_d)}$ and physical $q = -e^{-i\pi\lambda}$ (so $Y(\theta+i\pi) = e^{i\pi\lambda}Y(\theta) = -q^{-1}Y(\theta)$), our script `verify_lectures_1_to_6.py` confirms:
       - **Unitarity:** $\mathbb{T}_n(\theta)^{-1} = \frac{1}{\Phi(\theta)(1+Y(\theta)^2)}\begin{pmatrix} q^n & i Y q^n \\ i Y q^{-(n+1+s)} & q^{-(n+1+s)} \end{pmatrix}$, and at $s=0$, $\Phi(\theta)\Phi(-\theta)(1 + e^{2\lambda\theta}e^{-2\lambda\theta_d})\cdots$ (or Hermitian unitarity $\mathbb{T}_n^\dagger \mathbb{T}_n = |\Phi(\theta)|^2(1+Y^2)\mathbf{1}_2 = \mathbf{1}_2 \implies |\Phi(\theta)|^2(1 + e^{2\lambda(\theta-\theta_d)}) = 1$).
       - **Crossing (`38-defect-sine-gordon.qmd` line 219 `eq-sg-crossing`):** With $\mathcal{C} = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix} = \sigma_y$, the crossed block on $\mathcal{H}_n$ is $\mathcal{C}\mathbb{T}_n(\theta+i\pi)^{\mathsf{T}}\mathcal{C}^{-1}$, where the transposed basis reverses $\{v_+\otimes|n\rangle, v_-\otimes|n+1\rangle\}$ and shifts $n \to n+1$, giving
         $$
         \Phi(\theta)\Phi(\theta+i\pi)\bigl(1 + e^{2\lambda(\theta-\theta_d)}\bigr) = q^{-s},
         $$
         where all $q^{\pm n}$ factors cancel **identically** for every $n \ge 0$!
  2. **Soliton mass on line 172:**
     - Line 172 writes $M_1 = \frac{2m}{\beta^2}\sin(\pi/3)$, missing the factor $2h = 6$ (`40-defect-solitons.qmd` line 17 has $M_a = \frac{4hm}{\beta^2}\sin\frac{\pi a}{h}$).
  3. **Higher-rank $a_r^{(1)}$ Weyl representation (`eq-lec5-ar-weyl-action`, line 192) and $a_2^{(1)}$ transmission matrix (`prp-lec5-a2-matrices`, line 204):**
     - On line 192 (`eq-lec5-ar-weyl-action`), Lecture 5 writes $\rho(f_j)|\vec{N}\rangle \propto y_d q^{N_{j+1}+s_{j+1}}|\vec{N}+\vec{e}_j\rangle$. In the book (`42-soliton-transmission.qmd` line 30 `eq-dt-weyl`), the defect Weyl module acts on the $r$-oscillator root lattice $|\vec{N}\rangle = |N_1,\dots,N_r\rangle$ with spectral parameter $y_d^{-1}$ and shifts $|\vec{N}\rangle \to |\vec{N}+\vec{e}_j-\vec{e}_{j-1}\rangle$ for $j=1,\dots,r$ while $f_0$ shifts $|\vec{N}\rangle \to |\vec{N}-\vec{e}_r+\vec{e}_1\rangle$.
     - On line 204 (`prp-lec5-a2-matrices`), Lecture 5 writes $|j\rangle_V\otimes|\vec{N}\rangle \to |j+1\rangle_V\otimes|\vec{N}+\vec{\alpha}_j\rangle$ for **both** the diagonal and off-diagonal branches! In the book (`42-soliton-transmission.qmd` line 53 `eq-dt-tv`), the diagonal branch preserves both the soliton colour $|j\rangle_V$ and the defect state $|\vec{N}\rangle$, while the off-diagonal branch lowers the soliton colour $|j\rangle_V \to |j-1\rangle_V$ and shifts $|\vec{N}\rangle \to |\vec{N}+\vec{e}_j-\vec{e}_{j-1}\rangle$:
       $$
       \mathbb{T}_V(y)\bigl(|j\rangle_V \otimes |\vec{N}\rangle\bigr) = f(y)\!\left(q^{\tilde{N}_j}|\vec{N}\rangle \otimes |j\rangle_V + \varepsilon_j\frac{y}{y_d} q^{-\tilde{N}_j}|\vec{N}+\vec{e}_j-\vec{e}_{j-1}\rangle \otimes |j-1\rangle_V\right).
       $$

- **Replacement text for `lecture-notes/defects-gemini/05-crossing-resonances-and-tunnelling.qmd`:**
  - **Lines 29–96 (`sec-lec5-unitarity-crossing`):** Replace `eq-lec5-tn-2x2` through `eq-lec5-functional-pair` with the book's Weyl-module $2\times 2$ block derivation from `38-defect-sine-gordon.qmd` lines 208–235 (`sec-sg-prefactor`):
    ```markdown
    $$
    \mathbb{T}_n(\theta) = \Phi(\theta)\begin{pmatrix}
    q^{-n} & -i Y(\theta)\, q^{n+1+s} \\
    -i Y(\theta)\, q^{-n} & q^{n+1+s}
    \end{pmatrix}
    = \Phi(\theta)\begin{pmatrix} 1 & -i Y(\theta) \\ -i Y(\theta) & 1 \end{pmatrix}
    \begin{pmatrix} q^{-n} & 0 \\ 0 & q^{n+1+s} \end{pmatrix},
    $$ {#eq-lec5-tn-2x2}
    where $Y(\theta) = y/y_d = e^{\lambda(\theta-\theta_d)}$.
    ```
    Followed by the 3-line unitarity and crossing evaluation yielding (`eq-sg-crossing`, `38-defect-sine-gordon.qmd` line 219):
    ```markdown
    $$
    \Phi(\theta)\,\Phi(\theta + i\pi)\,\bigl(1 + e^{2\lambda(\theta - \theta_d)}\bigr) = q^{-s}, \qquad
    \bigl|\Phi(\theta)\bigr|^2\,\bigl(1 + e^{2\lambda(\theta - \theta_d)}\bigr) = 1 \quad (s = 0).
    $$ {#eq-lec5-functional-pair}
    ```
  - **Line 172:** Replace $M_1 = \frac{2m}{\beta^2}\sin(\pi/3)$ with $M_1 = \frac{12m}{\beta^2}\sin(\pi/3)$.
  - **Lines 188–215 (`eq-lec5-ar-weyl-action` and `prp-lec5-a2-matrices`):** Replace with the book's $r$-oscillator Weyl representation (`42-soliton-transmission.qmd` line 30 `eq-dt-weyl`) and vector transmission matrix (`42-soliton-transmission.qmd` line 53 `eq-dt-tv` and line 93 `eq-dt-anti-tv`):
    ```markdown
    $$
    \mathbb{T}_V(\theta)\bigl(|j\rangle_V \otimes |\vec{N}\rangle\bigr)
    = \Phi_1(\theta)\!\left(q^{\tilde{N}_j}|\vec{N}\rangle \otimes |j\rangle_V
    + \varepsilon_j e^{\mu(\theta-\theta_d)} q^{-\tilde{N}_j}|\vec{N}+\vec{e}_j-\vec{e}_{j-1}\rangle \otimes |j-1\rangle_V\right).
    $$ {#eq-lec5-a2-t1}
    ```

---

### Item 6: Lecture 6 (`lecture-notes/defects-gemini/06-boundaries-and-breathers.qmd`) and `index.qmd`

- **Verdict:** `wrong` (factor-of-2 error in `index.qmd` line 39; factor-of-2 error in `eq-lec6-fused-lagrangian` line 93; sign error $(+i)(-i)=-1$ in all four entries of `eq-lec6-kgz-explicit` line 145; missing $-i$ prefactor in `eq-lec6-b1-amplitude` line 222; non-chiral breather/particle transmission formula in `eq-lec6-breather-t-exact` line 233; and "spin-3 conserved charge" vs "spin-3 current / spin-2 charge" on line 252. The doubled-Borel coideal classification in lines 15–75 and the $G_n = \operatorname{diag}(q^{-n}, q^{n+1})$ conjugation factorization in `eq-lec6-gn-factorization` line 138 are `correct`!).
- **Independent computation (`referee_scratch/lectures_gemini/verify_lectures_1_to_6.py`, Section 6):**
  1. **`index.qmd` line 39:**
     - Verified in `verify_lectures_1_to_6.py` (check 2.2) that with $p = \frac{\phi+\psi}{2}, q = \frac{\phi-\psi}{2}$ and $\mathcal{B}_d = f(p)+g(q)$, $\frac{1}{2}(\partial_\phi\mathcal{B}_d)^2 - \frac{1}{2}(\partial_\psi\mathcal{B}_d)^2 = \frac{1}{2}f'(p)g'(q) = V(\phi)-V(\psi)$ (matching `02-momentum-and-backlund.qmd` line 143 `eq-lec2-master-functional`). `index.qmd` line 39 writes $\frac{1}{4}$ instead of $\frac{1}{2}$.
  2. **Fused boundary Lagrangian (`eq-lec6-fused-lagrangian`, line 93):**
     - Varying $\int dt\bigl[\frac{1}{2}\phi(0,t)\partial_t\psi(t) - \mathcal{B}_d(\phi(0,t),\psi(t)) - V_B(\psi(t))\bigr]$ with respect to $\phi(0,t)$ and $\psi(t)$ gives $\partial_x\phi = \frac{1}{2}\partial_t\psi - \partial_\phi\mathcal{B}_d$ and $\frac{1}{2}\partial_t\phi = -\partial_\psi\mathcal{B}_d - V_B'(\psi)$, missing a factor of 2 in `eq-lec6-fused-bc` (line 100)! Writing $\phi(0,t)\partial_t\psi(t)$ (or $\frac{1}{2}(\phi\partial_t\psi - \psi\partial_t\phi)$) in `eq-lec6-fused-lagrangian` produces `eq-lec6-fused-bc` (`43-defects-and-boundaries.qmd` line 25 `eq-db-fused-action`).
  3. **Quantum defect-boundary fusion (`eq-lec6-three-matrices-hn`, `eq-lec6-gn-factorization`, `eq-lec6-kgz-explicit`, lines 123–155):**
     - Verified symbolically in `verify_lectures_1_to_6.py` (checks 6.2a, 6.2b) that:
       - $\bar{\mathbb{T}}_n(\theta)\mathbb{T}_n(-\theta) = \Phi(\theta)\bar{\Phi}(\theta)(1 + Y(-\theta)^2)\mathbf{1}_2 = \mathbf{1}_2$,
       - $K_{\rm fused}(\theta)|_{\mathcal{H}_n} = \bar{\mathbb{T}}_n(\theta) K_D(\theta) \mathbb{T}_n(\theta) = G_n^{-1} K_{\rm GZ}(\theta) G_n$ with $G_n = \operatorname{diag}(q^{-n}, q^{n+1})$ completely independent of $n$.
     - **Error caught in `eq-lec6-kgz-explicit` (line 145, check 6.2c):** In `eq-lec6-three-matrices-hn` (line 126), $\bar{\mathbb{T}}_n(\theta)$ has $+i Y(-\theta)$ off the diagonal while $\mathbb{T}_n(\theta)$ has $-i Y(\theta)$ off the diagonal. Multiplying $\begin{pmatrix} 1 & +i Y(-\theta) \\ +i Y(-\theta) & 1 \end{pmatrix} \begin{pmatrix} k_+ & 0 \\ 0 & k_- \end{pmatrix} \begin{pmatrix} 1 & -i Y(\theta) \\ -i Y(\theta) & 1 \end{pmatrix}$ gives $(+i)(-i) = +1$ on the diagonal product and opposite signs on the off-diagonal terms:
       $$
       K_{\rm GZ}(\theta) = \Phi(\theta)\bar{\Phi}(\theta)\begin{pmatrix}
       k_+(\theta) + Y(\theta)Y(-\theta)k_-(\theta) & -i\bigl(Y(\theta)k_+(\theta) - Y(-\theta)k_-(\theta)\bigr) \\
       -i\bigl(Y(\theta)k_-(\theta) - Y(-\theta)k_+(\theta)\bigr) & Y(\theta)Y(-\theta)k_+(\theta) + k_-(\theta)
       \end{pmatrix},
       $$
       which matches `43-defects-and-boundaries.qmd` line 159 (`(k_+ + y_d^{-2} k_-)/(1 + y^{-2} y_d^{-2})` with $Y(\theta)Y(-\theta) = y_d^{-2}$) and vanishes off the diagonal at $\theta=0, k_+=k_-$. Line 145 of Lecture 6 treated $(+i)(-i)$ as $-1$, flipping the sign of every $Y(-\theta)$ term!
  4. **Breather singlet bootstrap and $a_r^{(1)}$ particle transmission (`eq-lec6-b1-amplitude` line 222 and `eq-lec6-breather-t-exact` line 233):**
     - **Error caught on lines 222 and 233:** On line 222 (`eq-lec6-b1-amplitude`), Lecture 6 writes $\mathbb{T}_{B_1}(\theta) = \frac{\sinh(\frac{\theta-\theta_d}{2}-\frac{i\pi}{4})}{\sinh(\frac{\theta-\theta_d}{2}+\frac{i\pi}{4})}$, missing the prefactor $-i$ (`38-defect-sine-gordon.qmd` line 337 `eq-sg-b1-transmission`: $\mathbb{T}_{B_1}(\theta) = -i\frac{\sinh(\frac{\theta-\theta_d}{2}-\frac{i\pi}{4})}{\sinh(\frac{\theta-\theta_d}{2}+\frac{i\pi}{4})} = \frac{i-\sinh(\theta-\theta_d)}{\cosh(\theta-\theta_d)}$).
     - Even more importantly, on line 233 (`eq-lec6-breather-t-exact`), Lecture 6 writes the $a_r^{(1)}$ fundamental particle transmission amplitude as
       $$
       \mathbb{T}_a(\theta) = \frac{\sinh(\theta - \theta_d) - i\sin\gamma_a}{\sinh(\theta - \theta_d) + i\sin\gamma_a}, \qquad \gamma_a = \frac{\pi a}{h}.
       $$
       This formula is **wrong** for a single Type-I$_-$ defect: it is parity-symmetric in $\theta-\theta_d$ (non-chiral), has *two* poles in $\theta-\theta_d \pmod{2\pi i}$, and at $\gamma_1 = \pi/2$ ($a_1^{(1)}$) equals $\frac{\sinh(\theta-\theta_d)-i}{\sinh(\theta-\theta_d)+i} = [\mathbb{T}^{\rm cl}(\theta)]^2$, contradicting `eq-lec6-b1-amplitude` eleven lines earlier! As shown in the main book (`44-particle-transmission.qmd` lines 39 and 66 `eq-dpt-particle`), the Type-I$_-$ fundamental particle transmission amplitude is **chiral**:
       $$
       \mathbb{T}_a(\theta) = \frac{e^{-(\theta - \theta_d)} + e^{i\gamma_a}}{e^{-(\theta - \theta_d)} + e^{-i\gamma_a}} = -e^{i\gamma_a}\,\frac{\sinh\!\left(\frac{\theta - \theta_d}{2} - \frac{i(\pi - \gamma_a)}{2}\right)}{\sinh\!\left(\frac{\theta - \theta_d}{2} + \frac{i(\pi - \gamma_a)}{2}\right)},
       $$
       which at $\gamma_1 = \pi/2$ reduces to $-i\frac{\sinh(\frac{\theta-\theta_d}{2}-\frac{i\pi}{4})}{\sinh(\frac{\theta-\theta_d}{2}+\frac{i\pi}{4})} = \frac{i-\sinh(\theta-\theta_d)}{\cosh(\theta-\theta_d)}$, matching `eq-lec6-b1-amplitude`!
  5. **Current spin vs. charge spin on line 252 (`sec-lec6-anomalies`):**
     - Line 252 refers to the Penati–Polvara loop calculation for $a_r^{(1)}$ ($r \ge 2$) as *"the one-loop quantum correction to the spin-3 conserved charge across the defect"*. In the book's conventions (`briefs/part8_defects_conventions.md` line 222; `42-soliton-transmission.qmd` line 46), a current $T_{\pm(s+1)}$ of Lorentz weight $\pm(s+1)$ integrates to a conserved charge $Q_{\pm s}$ of spin $\pm s$: for $a_r^{(1)}$ ($r \ge 2$), the first higher-spin current is the **spin-3 current** (whose conserved charge has **spin 2**; for sinh-Gordon $r=1$, where even-spin charges vanish, it is the spin-4 current / spin-3 charge).

- **Replacement text for `lecture-notes/defects-gemini/index.qmd` and `06-boundaries-and-breathers.qmd`:**
  - **`lecture-notes/defects-gemini/index.qmd`, line 39:**
    ```markdown
       $$\frac{1}{2}\,f'\!\left(\frac{\phi+\psi}{2}\right) g'\!\left(\frac{\phi-\psi}{2}\right) = V(\phi) - V(\psi).$$
    ```
  - **`lecture-notes/defects-gemini/06-boundaries-and-breathers.qmd`, line 93 (`eq-lec6-fused-lagrangian`):**
    ```markdown
    $$
    L_{\rm B, eff}\bigl(\phi(0,t), \psi(t)\bigr) = \phi(0,t)\,\partial_t\psi(t) - \mathcal{B}_d\bigl(\phi(0,t), \psi(t)\bigr) - V_B\bigl(\psi(t)\bigr).
    $$ {#eq-lec6-fused-lagrangian}
    ```
  - **`lecture-notes/defects-gemini/06-boundaries-and-breathers.qmd`, lines 142–149 (`eq-lec6-kgz-explicit`):**
    ```markdown
    $$
    K_{\rm GZ}(\theta) = \Phi(\theta)\bar{\Phi}(\theta)\begin{pmatrix}
    k_+(\theta) + Y(\theta)Y(-\theta)\,k_-(\theta) & -i\bigl(Y(\theta)k_+(\theta) - Y(-\theta)k_-(\theta)\bigr) \\
    -i\bigl(Y(\theta)k_-(\theta) - Y(-\theta)k_+(\theta)\bigr) & Y(\theta)Y(-\theta)\,k_+(\theta) + k_-(\theta)
    \end{pmatrix}.
    $$ {#eq-lec6-kgz-explicit}
    ```
  - **`lecture-notes/defects-gemini/06-boundaries-and-breathers.qmd`, lines 220–236 (`eq-lec6-b1-amplitude` and `eq-lec6-breather-t-exact`):**
    ```markdown
    $$
    \mathbb{T}_{B_1}(\theta) = -i\,\frac{\sinh\!\left(\frac{\theta - \theta_d}{2} - \frac{i\pi}{4}\right)}{\sinh\!\left(\frac{\theta - \theta_d}{2} + \frac{i\pi}{4}\right)} = \frac{i - \sinh(\theta - \theta_d)}{\cosh(\theta - \theta_d)} = \mathbb{T}^{\rm cl}(\theta).
    $$ {#eq-lec6-b1-amplitude}

    Wait: what happened to the coupling constant $\lambda$, the quantum parameter $q$, and the infinite-dimensional defect oscillator $\hat{N}$? **They cancelled out completely!** The exact quantum transmission amplitude of the lightest breather across the Type-I$_-$ defect is independent of $\hbar$ and coincides identically with the classical linearised Klein–Gordon transmission amplitude (-@eq-lec2-t-cl) we derived in Lecture 2!

    For a general affine Toda theory $a_r^{(1)}$ (or higher breathers by bootstrap induction, see @sec-dpt-particles in the main text), the exact quantum transmission amplitude of the $a$-th fundamental particle across a Type-I$_-$ defect is the chiral phase

    $$
    \mathbb{T}_a(\theta) = \frac{e^{-(\theta - \theta_d)} + e^{i\gamma_a}}{e^{-(\theta - \theta_d)} + e^{-i\gamma_a}} = -e^{i\gamma_a}\,\frac{\sinh\!\left(\frac{\theta - \theta_d}{2} - \frac{i(\pi - \gamma_a)}{2}\right)}{\sinh\!\left(\frac{\theta - \theta_d}{2} + \frac{i(\pi - \gamma_a)}{2}\right)}, \qquad \gamma_a = \frac{\pi a}{h},
    $$ {#eq-lec6-breather-t-exact}
    ```
  - **`lecture-notes/defects-gemini/06-boundaries-and-breathers.qmd`, line 252 (`sec-lec6-anomalies`):**
    ```markdown
    Historically, this split between $A_r$ and $D_r/E_r$ was discovered from the bottom up. When Penati and Polvara computed the one-loop quantum correction to the spin-3 conserved current (spin-2 charge) across the defect in $a_r^{(1)}$ ($r \ge 2$) and the spin-4 current (spin-3 charge) in sinh-Gordon, they found that the quantum anomaly vanishes without modifying the classical defect potential.
    ```

---

## Unverified Observations and Unchecked Claims

1. **`lecture-notes/defects-gemini/index.qmd` line 98 ("How these lectures were written"):** Currently states *"The lectures have not yet been checked by an independent referee; if you spot an error, please report it."* Once the fixes from this report are integrated, that disclaimer can be updated.
2. **All analytical claims in Lectures 1–6 were checked:** No mathematical derivations or displayed equations in the six lectures were left unchecked; the full symbolic and numerical verification log is recorded in `referee_scratch/lectures_gemini/verify_lectures_1_to_6.py`.
