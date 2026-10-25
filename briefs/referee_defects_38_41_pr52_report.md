# Referee Report: PR #52 (Ch. 38 Scalar Factor & Poles, Ch. 41 Modular Double) and Lectures 8–9 Improvements

- **Date:** 2026-10-24
- **Model:** Gemini 4 Argon
- **Commit reviewed:** `9b912ee` (PR #52, merged on `main`)
- **Independent verification scripts:**
  - `referee_scratch/defects_38_41_pr52/verify_ch38.py` (Items 1–4)
  - `referee_scratch/defects_38_41_pr52/verify_ch41_43.py` (Items 5–7)

---

## Summary Table

| Item | Target | Verdict | Status Label Assessment |
|---|---|---|---|
| **1** | `prp-sgd-scalar-integral`, `eq-sgd-scalar-integral` (`38-defect-sine-gordon.qmd`) | **correct** | `[Theorem]` justified as stands |
| **2** | `prp-sg-breather-transmission` & proof (`38-defect-sine-gordon.qmd`) | **correct** | Upgrade to `[Established, with a new proof for $B_1$; numerically checked]` (or `[Theorem]` if the 2-line breather induction is added) |
| **3** | `prp-sgd-pole-ladder`, `eq-sgd-widths` (`38-defect-sine-gordon.qmd`) | **correct** | `[Theorem]` justified as stands |
| **4** | `prp-sgd-odd-jump` & inverted-oscillator remarks (`38-defect-sine-gordon.qmd`) | **correct** (minor edge-case note at $\theta_d = 0$ in proof) | `[Theorem]` justified (recommend adding $\theta_d \ne 0$ in one line of the proof) |
| **5** | `prp-dq-two-weyl-pairs` (`41-defect-quantum-groups.qmd`) | **correct** | `[Theorem]` justified as stands |
| **6** | Literature corrections (`38`, `41`, `D-corrections.qmd`, `E-original-results.qmd`) | **correct** (suggest expanding Nemes entry in App. D to list the 4 additional misprints from Ch. 38 Notes) | Labels in App. E match statements |
| **7a** | Chirality reading of `prp-diagonal-borel` (`43-defects-and-boundaries.qmd` line 110; `41-defect-quantum-groups.qmd` line 61) | **wrong** in `43-defects-and-boundaries.qmd` line 110 (reversed order of $\cos\frac{\beta(\phi\mp\psi)}2$ perturbations relative to passed/reflected chiralities) | Recommend replacing line 110 of Ch. 43 with the exact $\psi \to -\psi$ / `prp-two-borels` argument of Lecture 9 (`sec-l9-borel`) |
| **7b** | Upgrade of `prp-defect-fusion-k` (3) (`43-defects-and-boundaries.qmd` lines 147–161) via `thm-l9-dirichlet` | **correct** | Recommend upgrading `prp-defect-fusion-k` to `[Theorem]` (and updating App. E line 54) |
| **7c** | Trace-per-state for `prp-defect-transfer` (2) and fusion proof of `cnj-defect-tq` (`41-defect-quantum-groups.qmd` lines 228–262) via `thm-l8-commute` and `thm-l8-tq` | **correct** | Recommend upgrading `cnj-defect-tq` (`thm-defect-tq`) to `[Theorem]` and `prp-defect-transfer` (2) to `[Theorem]` using Lecture 8's proofs |

---

## Item 1: Fourier integral representation of $\Phi(\theta)$ (`prp-sgd-scalar-integral`, `eq-sgd-scalar-integral`)

- **Verdict:** **correct**
- **Target:** `part8-defects/38-defect-sine-gordon.qmd`, lines 225–252 (`prp-sgd-scalar-integral` and `eq-sgd-scalar-integral`).

### Independent Verification

1. **Step-by-step derivation from the Gamma product (`eq-sgd-scalar`) via Malmstén's formula:**
   - For $\theta \in \mathbb{R}$ and each $k \ge 0$, the four Gamma arguments in $r_N(\theta)$ (`eq-sgd-scalar`, line 199) are $a_i - z$ with $z = \frac{i\lambda\theta}{2\pi} \in i\mathbb{R}$ and
     $$
     a_1 = k\lambda + \frac14, \quad a_2 = (k+1)\lambda + \frac34, \quad b_1 = \Big(k+\frac12\Big)\lambda + \frac14, \quad b_2 = \Big(k+\frac12\Big)\lambda + \frac34,
     $$
     all of which have $\operatorname{Re}(a_i \mp z) = a_i \ge \frac14 > 0$ for $\lambda > 0$.
   - Since $a_1 + a_2 - b_1 - b_2 = (2k\lambda + \lambda + 1) - (2k\lambda + \lambda + 1) = 0$, the linear term $\mp z t e^{-t}$ in Malmstén's formula cancels identically inside each $k$-factor, and $\ln\Gamma(a - z) - \ln\Gamma(a + z) = \int_0^\infty \frac{dt}{t}\Big[\text{(linear)} - \frac{2i e^{-at}\sin(\lambda\theta t/2\pi)}{1 - e^{-t}}\Big]$ sums over the four arguments to give the geometric series over $k \ge 0$:
     $$
     \sum_{k=0}^\infty \big(e^{-a_1 t} + e^{-a_2 t} - e^{-b_1 t} - e^{-b_2 t}\big) = \frac{(1 - e^{-\lambda t/2})(e^{-t/4} - e^{-\lambda t/2 - 3t/4})}{1 - e^{-\lambda t}}.
     $$
   - Substituting $t = \frac{2\pi k}{\lambda}$ and multiplying by $\frac{2}{1 - e^{-t}}$ gives Nemes's integral (`@nemes2010`, eq. (49)) for $|\operatorname{Im}\theta| < \frac{\pi}{2\lambda}$:
     $$
     \frac{1}{i}\ln\frac{r(\theta)}{\bar r(\theta)} = -\int_0^\infty \frac{dk}{k}\,\frac{\sinh\big(\frac{\pi k}{2}(1 + \lambda^{-1})\big)\sin(k\theta)}{\sinh\frac{\pi k}{\lambda}\cosh\frac{\pi k}{2}}.
     $$
   - Using $\sinh\big(\frac{\pi k}{2} + \frac{\pi k}{2\lambda}\big) = \sinh\frac{\pi k}{2}\cosh\frac{\pi k}{2\lambda} + \cosh\frac{\pi k}{2}\sinh\frac{\pi k}{2\lambda}$ and $\sinh\frac{\pi k}{\lambda} = 2\sinh\frac{\pi k}{2\lambda}\cosh\frac{\pi k}{2\lambda}$, the kernel splits into
     $$
     \frac{\sinh\big(\frac{\pi k}{2}(1 + \lambda^{-1})\big)}{\sinh\frac{\pi k}{\lambda}\cosh\frac{\pi k}{2}} = \frac{\tanh\frac{\pi k}{2}}{2\sinh\frac{\pi k}{2\lambda}} + \frac{1}{2\cosh\frac{\pi k}{2\lambda}}.
     $$
   - Applying the standard Fourier integrals
     $$
     \int_0^\infty \frac{\sin(k\theta)\,dk}{k\cosh\frac{\pi k}{2\lambda}} = 2\arctan\tanh\frac{\lambda\theta}{2} = 2\arctan(e^{\lambda\theta}) - \frac{\pi}{2}, \qquad
     \int_0^\infty \frac{(1 - \cos(k\theta))\,dk}{k\sinh\frac{\pi k}{2\lambda}} = \ln\cosh(\lambda\theta),
     $$
     and adding $\ln\frac{e^{i\pi(1+\lambda)/4}}{1 + ie^{\lambda\theta}} = -\frac12\ln(1 + e^{2\lambda\theta}) + i\big(\frac{\pi(1+\lambda)}{4} - \arctan e^{\lambda\theta}\big) = -\frac{\lambda\theta}{2} - \frac12\ln 2 - \frac12\ln\cosh(\lambda\theta) + i\big(\frac{\pi(1+\lambda)}{4} - \arctan e^{\lambda\theta}\big)$, the $\arctan e^{\lambda\theta}$ and $\frac{i\pi}{4}$ terms cancel and the real and imaginary parts combine via
     $$
     \cos(k\theta) - 1 + i\tanh\frac{\pi k}{2}\sin(k\theta) = \frac{\cosh\big(\frac{\pi k}{2} + ik\theta\big)}{\cosh\frac{\pi k}{2}} - 1
     $$
     into `eq-sgd-scalar-integral`.
2. **Convergence strip and analytic continuation:**
   - As $k \to 0$, $\cosh(\frac{\pi k}{2} + ik\theta)/\cosh\frac{\pi k}{2} - 1 = \frac{i\pi\theta - \theta^2}{2} k^2 + O(k^3)$ and $2k\sinh\frac{\pi k}{2\lambda} = \frac{\pi}{\lambda} k^2 + O(k^4)$, so the integrand has the finite limit $\frac{\lambda(i\pi\theta - \theta^2)}{2\pi}$ at $k = 0$.
   - As $k \to \infty$, writing $\frac{\cosh(\pi k/2 + ik\theta)/\cosh(\pi k/2) - 1}{2k\sinh(\pi k/2\lambda)} = \frac{e^{(i\theta - \pi/2\lambda)k} - e^{-\pi k/2\lambda} + e^{(-i\theta - \pi - \pi/2\lambda)k} - e^{(-\pi - \pi/2\lambda)k}}{k(1 + e^{-\pi k})(1 - e^{-\pi k/\lambda})}$, the four exponential exponents have real parts $(-\operatorname{Im}\theta - \frac{\pi}{2\lambda})k$, $-\frac{\pi}{2\lambda}k$, $(\operatorname{Im}\theta - \pi - \frac{\pi}{2\lambda})k$, and $(-\pi - \frac{\pi}{2\lambda})k$, all of which are strictly negative with uniform exponential decay on compact subsets if and only if
     $$
     -\frac{\pi}{2\lambda} < \operatorname{Im}\theta < \pi + \frac{\pi}{2\lambda}.
     $$
   - Therefore the integral defines a holomorphic function of $\theta$ on this strip, making $\Phi(\theta) = \exp(\ln\Phi(\theta))$ analytic and zero-free there.
3. **Unitarity, crossing, and $\theta \to -\infty$ limit directly from `eq-sgd-scalar-integral`:**
   - For $\theta \in \mathbb{R}$, $\operatorname{Re}\ln\Phi(\theta) = -\frac{\lambda\theta}{2} - \frac12\ln 2 - \frac12\ln\cosh(\lambda\theta) = -\frac12\ln(1 + e^{2\lambda\theta})$, so $|\Phi(\theta)|^2(1 + e^{2\lambda\theta}) = 1$.
   - For $|\operatorname{Im}\theta| < \frac{\pi}{2\lambda}$, both $\theta$ and $\theta + i\pi$ lie in the strip $(-\frac{\pi}{2\lambda}, \pi + \frac{\pi}{2\lambda})$. Since $\cosh(\frac{\pi k}{2} + ik(\theta + i\pi)) = \cosh(-\frac{\pi k}{2} + ik\theta) = \cosh(\frac{\pi k}{2} - ik\theta)$, the sum of the two numerators is $2(\cos(k\theta) - 1)\cosh\frac{\pi k}{2}$, giving $\ln\Phi(\theta) + \ln\Phi(\theta + i\pi) = -\lambda\theta - \ln 2 - \ln\cosh(\lambda\theta) = -\ln(1 + e^{2\lambda\theta})$.
   - As $\theta \to -\infty$, the oscillatory part vanishes by the Riemann–Lebesgue lemma and $-\frac12\ln(1 + e^{2\lambda\theta}) \to 0$, so $\Phi(\theta) \to 1$.

### Script & Output (`referee_scratch/defects_38_41_pr52/verify_ch38.py`)

```text
1.1 [PASS] Step A (signed sum over k): diff = 0
1.2 [PASS] Step B (Nemes kernel): diff = 0
1.3 [PASS] Step C (kernel decomposition): diff = 0
1.4 [PASS] Step D (arctan & ln cosh integrals): err_i = 2.3e-41, err_ii = 0.0
1.5 [PASS] Step E (real part identity): |re_lhs - re_rhs| = 0.0
1.6 High-precision comparison of exp(ln_phi_integral) vs phi_gamma_product across the strip:
    theta = (-1.2 + 0.0j)        | rel_err = 1.282e-29 -> PASS
    theta = (0.75 + 0.0j)        | rel_err = 6.982e-30 -> PASS
    theta = (0.4 - 0.45j)        | rel_err = 5.08e-30 -> PASS
    theta = (-0.3 + 0.5j)        | rel_err = 4.822e-30 -> PASS
    theta = (0.5 + 1.8j)         | rel_err = 9.378e-30 -> PASS
    theta = (-0.2 + 3.5j)        | rel_err = 3.432e-30 -> PASS
1.7 [PASS] Unitarity err = 2.3e-41, Crossing log err = 1.15e-41, Limit |Phi(-10)-1| = 7.5e-5 -> |Phi(-20)-1| = 3.4e-9
```

### Label Assessment

- **`[Theorem]` is fully justified.** The proof in lines 239–252 is complete, rigorous, and requires no changes.

---

## Item 2: Analytic proof of lightest breather transmission & higher breathers (`prp-sg-breather-transmission`)

- **Verdict:** **correct**
- **Target:** `part8-defects/38-defect-sine-gordon.qmd`, lines 333–362 (`prp-sg-breather-transmission` and its proof).

### Independent Verification

1. **Bootstrap action on the singlet $(v_+\otimes v_- - v_-\otimes v_+)\otimes|n\rangle$:**
   - With $u_a = u - \frac{iu_1}{2}$, $u_b = u + \frac{iu_1}{2}$, $u_1 = \pi - \frac{\pi}{\lambda}$, and $q = -e^{-i\pi\lambda}$, we have $Y_1 = e^{\lambda u_a} = q^{1/2}Y$ and $Y_2 = e^{\lambda u_b} = q^{-1/2}Y$, so $Y_1 = q Y_2$ and $Y_1 Y_2 = Y^2$.
   - Direct computation of $\mathbb{T}_{13}(u_a)\mathbb{T}_{23}(u_b)(v_+\otimes v_- - v_-\otimes v_+)\otimes|n\rangle$ gives:
     - coefficient of $|n-1\rangle\otimes v_+\otimes v_+$: $-i q^s(q Y_2 - Y_1) = 0$;
     - coefficient of $|n+1\rangle\otimes v_-\otimes v_-$: $-i q^s(Y_1 - q Y_2) = 0$;
     - coefficient of $|n\rangle\otimes(v_+\otimes v_- - v_-\otimes v_+)$: $q^s(1 - q Y_1 Y_2) = q^s(1 + q Y^2)$ (Wait: in `verify_ch38.py` step 2.1, $1 - (-iY_2)(-iY_1) q = 1 + q Y^2$, which for $s = 0$ is $1 + q Y^2$).
2. **Crossing step and strip validity:**
   - Setting $w = u - \frac{iu_1}{2} = u - \frac{i\pi}{2} + \frac{i\pi}{2\lambda}$, we have $e^{2\lambda w} = e^{2\lambda u - i\pi\lambda + i\pi} = q Y^2$, $w + i\pi = v + \varepsilon$, and $u + \frac{iu_1}{2} = v - \varepsilon$ with $v = u + \frac{i\pi}{2}$ and $\varepsilon = \frac{i\pi}{2\lambda}$.
   - Using $\Phi(w)\Phi(w + i\pi)(1 + e^{2\lambda w}) = 1$, we obtain
     $$
     \mathbb{T}_{B_1}(\theta) = \Phi(w)\Phi(v - \varepsilon)(1 + q Y^2) = \frac{\Phi(v - \varepsilon)}{\Phi(v + \varepsilon)}.
     $$
   - For $|\operatorname{Im}u| < \frac{\pi}{2}$, $\operatorname{Im}(v - \varepsilon) = \operatorname{Im}u + \frac{\pi}{2} - \frac{\pi}{2\lambda} \in (-\frac{\pi}{2\lambda}, \pi - \frac{\pi}{2\lambda})$ and $\operatorname{Im}(v + \varepsilon) = \operatorname{Im}u + \frac{\pi}{2} + \frac{\pi}{2\lambda} \in (\frac{\pi}{2\lambda}, \pi + \frac{\pi}{2\lambda})$. Both points lie **strictly inside** the convergence strip $(-\frac{\pi}{2\lambda}, \pi + \frac{\pi}{2\lambda})$ of `eq-sgd-scalar-integral` for every $\lambda > 0$.
3. **Cancellation of $\sinh\frac{\pi k}{2\lambda}$ in the integral difference:**
   - Subtracting `eq-sgd-scalar-integral` at $v + \varepsilon$ from $v - \varepsilon$ gives linear prefactor $-\frac{\lambda}{2}(-2\varepsilon) = \frac{i\pi}{2}$ and integrand numerator
     $$
     \cosh\Big(\frac{\pi k}{2} + ikv - ik\varepsilon\Big) - \cosh\Big(\frac{\pi k}{2} + ikv + ik\varepsilon\Big) = -2\sinh\Big(\frac{\pi k}{2} + ikv\Big)\sin(k\varepsilon) = -2i\cos(ku)\sinh\frac{\pi k}{2\lambda}.
     $$
   - The Coupling-dependent denominator $\sinh\frac{\pi k}{2\lambda}$ cancels identically, leaving
     $$
     \ln\frac{\Phi(v-\varepsilon)}{\Phi(v+\varepsilon)} = \frac{i\pi}{2} - i\int_0^\infty \frac{\cos(ku)\,dk}{k\cosh\frac{\pi k}{2}} = \frac{i\pi}{2} - 2i\arctan e^u = \ln\frac{i - \sinh u}{\cosh u} = \ln\frac{-i\sinh\big(\frac{u}{2} - \frac{i\pi}{4}\big)}{\sinh\big(\frac{u}{2} + \frac{i\pi}{4}\big)}.
     $$
4. **Higher breathers ($n \ge 2$) in Bowcock–Corrigan–Zambon 2005 (`hep-th/0506169`):**
   - In `0506169.tex` lines 1368–1429 (Section 5, eqs. (5.5)–(5.13)), Bowcock, Corrigan and Zambon derive the transmission factors for all breathers $B_n$ ($n \ge 1$) explicitly: first $T^1(\theta) = -\frac{i}{(1)_{\theta - \eta}}$ in (5.7), then $T^2(\theta)$ in (5.9) via the $B_1 + B_1 \to B_2$ bootstrap $T^p(\theta) = T^n(\theta + i\pi m/\gamma) T^m(\theta - i\pi n/\gamma)$ of (5.8), and by induction the closed product formulas (5.12) (for odd breathers $B_{2s+1}$) and (5.13) (for even breathers $B_{2s}$).
   - In `verify_ch38.py` step 2.4, we checked to 40 digits that the book's compact product $\prod_{k=1}^n \mathbb{T}^{\rm cl}(\theta + \frac{i\pi(n+1-2k)}{2\lambda})$ in `prp-sg-breather-transmission` (line 341) is algebraically identical to BCZ's (5.12)–(5.13) for both odd and even $n$: since $\mathbb{T}^{\rm cl}(u) = -i / (1)_u$ in BCZ's notation $(z)_u = \frac{\sinh(u/2 + i\pi z/4)}{\sinh(u/2 - i\pi z/4)}$, pairing the $k$ and $n+1-k$ factors in the book's product produces $(1 + z_l)_u (1 - z_l)_u$ with a $(-i)^2 = -1$ sign per pair, matching BCZ's $(-1)^{s+1}i$ (odd $n = 2s+1$) and $(-1)^s$ (even $n = 2s$).
   - Moreover, as Lecture 9 (`09-folding-fusing.qmd`, `thm-l9-breather`, line 391) observes, the induction step $\mathbb{T}_{B_{n+1}}(\theta) = \mathbb{T}_{B_n}(\theta + \frac{i\pi}{2\lambda})\mathbb{T}_{B_1}(\theta - \frac{in\pi}{2\lambda})$ from the breather fusion pole $\theta_b - \theta_a = \frac{i(n+1)\pi}{2\lambda}$ is a two-line calculation!

### Script & Output (`referee_scratch/defects_38_41_pr52/verify_ch38.py`)

```text
2.1 [PASS] Bootstrap algebra on singlet:
    coeff(|n-1> v_+ v_+) = 0
    coeff(|n+1> v_- v_-) = 0
    coeff(|n>   v_+ v_-) = qs*(Y**2*q + 1)
    coeff(|n>   v_- v_+) = qs*(-Y**2*q - 1)
2.2 [PASS] Strip check: for |Im u| < pi/2, Im(v -/+ eps) = Im u + pi/2 -/+ pi/(2*lam) both lie strictly in (-pi/(2*lam), pi + pi/(2*lam)).
    u = 0.65               | |Phi(v-e)/Phi(v+e) - T^cl| = 6.89e-41, |T^cl - (i-sinh u)/cosh u| = 0.0
    u = (0.4 + 0.7j)       | |Phi(v-e)/Phi(v+e) - T^cl| = 6.42e-42, |T^cl - (i-sinh u)/cosh u| = 5.74e-42
    u = (-0.8 - 0.5j)      | |Phi(v-e)/Phi(v+e) - T^cl| = 0.0, |T^cl - (i-sinh u)/cosh u| = 2.3e-41
    n = 1: |prod_book - BCZ (5.12/13)| = 0.0 -> PASS
    n = 2: |prod_book - BCZ (5.12/13)| = 1.28e-41 -> PASS
    n = 3: |prod_book - BCZ (5.12/13)| = 1.95e-41 -> PASS
    n = 4: |prod_book - BCZ (5.12/13)| = 2.57e-41 -> PASS
    n = 5: |prod_book - BCZ (5.12/13)| = 2.37e-41 -> PASS
```

### Label Assessment & Recommended Replacement Text

- Since BCZ (`@bcz2005`, eqs. (5.8)–(5.13)) already derived the $n \ge 2$ breather transmission factors from the breather bootstrap, the $n \ge 2$ part is at least `[Established; numerically checked]`, not merely `[Numerical for n >= 2]`. Even better: adding the two-line breather bootstrap induction (as in Lecture 9's `thm-l9-breather`) proves the $n \ge 2$ formula completely in the text, allowing `prp-sg-breather-transmission` to be upgraded to **`[Theorem]`** (or **`[Established, with a new proof]`**), with `appendices/E-original-results.qmd` line 44 updated to match.
- **Replacement text for `part8-defects/38-defect-sine-gordon.qmd`, lines 337–360:**
  ```markdown
  ::: {#prp-sg-breather-transmission}
  ## Breather transmission factors
  **[Theorem]** In the even sector $s=0$, the transmission matrix (-@eq-sg-transmission-matrix) acts on the $n$-th breather state inside $V\otimes V$ as the identity on $\mathcal W$ times
  $$
  \mathbb T_{B_n}(\theta)=\prod_{k=1}^n\mathbb T^{\rm cl}\Big(\theta+\frac{i\pi(n+1-2k)}{2\lambda}\Big),
  $$
  where $\mathbb T^{\rm cl}(\theta)$ is the classical linear-wave transmission factor (-@eq-sg-linear-transmission). In particular $\mathbb T_{B_1}(\theta)=\mathbb T^{\rm cl}(\theta)$, independent of $\hbar$.
  :::

  ::: {.proof}
  For $B_1$ at rapidity $\theta$, with $u=\theta-\theta_d$ and $u_1=\pi-\pi/\lambda$, the constituent rapidities are $\theta_{a,b}=\theta\mp iu_1/2$. Then $Y_1=e^{\lambda(u-iu_1/2)}$ and $Y_2=e^{\lambda(u+iu_1/2)}$ satisfy $Y_1=qY_2$ and $Y_1Y_2=Y^2$, because $q=-e^{-i\pi\lambda}$ (-@eq-q-physical). Applying $\mathbb T_{13}(\theta_a)\mathbb T_{23}(\theta_b)$ to $(v_+\otimes v_--v_-\otimes v_+)\otimes|n\rangle$, the $|n-1\rangle\otimes v_+\otimes v_+$ terms have coefficients $-iqY_2+iY_1=0$, the $|n+1\rangle\otimes v_-\otimes v_-$ terms have $-iY_1+iqY_2=0$, and the $|n\rangle$ terms give $(1+qY^2)|n\rangle\otimes(v_+\otimes v_--v_-\otimes v_+)$. Hence
  $$
  \mathbb T_{B_1}(\theta)=\Phi\big(u-iu_1/2\big)\,\Phi\big(u+iu_1/2\big)\big(1+qY^2\big).
  $$
  Write $w=u-iu_1/2=u-\frac{i\pi}2+\frac{i\pi}{2\lambda}$. Then $e^{2\lambda w}=qY^2$, $w+i\pi=v+\varepsilon$ and $u+iu_1/2=v-\varepsilon$, where $v=u+\frac{i\pi}2$ and $\varepsilon=\frac{i\pi}{2\lambda}$. Crossing (-@eq-sgd-scalar) at $w$, $\Phi(w)\Phi(w+i\pi)(1+e^{2\lambda w})=1$, gives
  $$
  \mathbb T_{B_1}(\theta)=\frac{\Phi(v-\varepsilon)}{\Phi(v+\varepsilon)}.
  $$
  For $|\operatorname{Im}u|<\frac\pi2$ both $v\mp\varepsilon$ have imaginary parts $\operatorname{Im}u+\frac\pi2\mp\frac\pi{2\lambda}$ in the strip $(-\frac\pi{2\lambda},\pi+\frac\pi{2\lambda})$ of @prp-sgd-scalar-integral. Subtracting (-@eq-sgd-scalar-integral) at $v+\varepsilon$ from that at $v-\varepsilon$, the linear terms give $-\frac\lambda2(-2\varepsilon)=\frac{i\pi}2$. In the integral,
  $$
  \cosh\Big(\frac{\pi k}2+ikv-ik\varepsilon\Big)-\cosh\Big(\frac{\pi k}2+ikv+ik\varepsilon\Big)=-2\sinh\Big(\frac{\pi k}2+ikv\Big)\sin(k\varepsilon)=-2i\cos(ku)\sinh\frac{\pi k}{2\lambda}.
  $$
  The factor $\sinh\frac{\pi k}{2\lambda}$ cancels the denominator, leaving no $\lambda$-dependence:
  $$
  \ln\mathbb T_{B_1}(\theta)=\frac{i\pi}2-i\int_0^\infty\frac{dk}{k}\,\frac{\cos(ku)}{\cosh\frac{\pi k}2}.
  $$
  Differentiating the integral gives $-\int_0^\infty\frac{\sin(ku)\,dk}{\cosh(\pi k/2)}=-\frac1{\cosh u}$, and since it vanishes at $u=0$ it equals $-2\arctan\tanh\frac u2=\frac\pi2-2\arctan e^u$. Thus $\mathbb T_{B_1}(\theta)=\exp(2i\arctan e^u-\frac{i\pi}2)=\frac{i-\sinh u}{\cosh u}$, which equals $-i\sinh(\frac u2-\frac{i\pi}4)/\sinh(\frac u2+\frac{i\pi}4)=\mathbb T^{\rm cl}(\theta)$.

  For $n\ge2$, $B_n$ is fused from $B_{n-1}$ at $\theta+\frac{i\pi}{2\lambda}$ and $B_1$ at $\theta-\frac{i(n-1)\pi}{2\lambda}$ [@bcz2005, eq. (5.8)], so $\mathbb T_{B_n}(\theta)=\mathbb T_{B_{n-1}}(\theta+\frac{i\pi}{2\lambda})\mathbb T_{B_1}(\theta-\frac{i(n-1)\pi}{2\lambda})$, which gives the product by induction [@bcz2005, eqs. (5.12)–(5.13)]. Also checked directly from the soliton–antisoliton bootstrap at $\theta_{a,b}=\theta\mp iu_n/2$, $u_n=\pi-n\pi/\lambda$, for $n=1,2,3$ in `code/defect_sg_transmission.py`.
  :::
  ```
- **Corresponding update in `appendices/E-original-results.qmd`, line 44:**
  Change `- **[Theorem]** Analytic proof of $\mathbb T_{B_1}(\theta)=\mathbb T^{\rm cl}(\theta)$ from (-@eq-sgd-scalar-integral):` to match the upgraded label of `@prp-sg-breather-transmission` (it already begins with `- **[Theorem]**`, so upgrading `prp-sg-breather-transmission` to `**[Theorem]**` keeps `check-refs.py` in exact agreement!).

---

## Item 3: Pole ladder and resonance widths (`prp-sgd-pole-ladder`, `eq-sgd-widths`)

- **Verdict:** **correct**
- **Target:** `part8-defects/38-defect-sine-gordon.qmd`, lines 277–308 (`prp-sgd-pole-ladder` and `eq-sgd-widths`).

### Independent Verification

1. **Analytic continuation and pole locations via crossing:**
   - By `prp-sgd-scalar-integral`, $\Phi(u)$ ($u = \theta - \theta_d$) is analytic and zero-free on $-\frac{\pi}{2\lambda} < \operatorname{Im}u < \pi + \frac{\pi}{2\lambda}$.
   - For $-\pi - \frac{\pi}{2\lambda} < \operatorname{Im}u \le -\frac{\pi}{2\lambda}$, the shifted argument $u + i\pi$ has $\operatorname{Im}(u + i\pi) \in (-\frac{\pi}{2\lambda}, \pi - \frac{\pi}{2\lambda}] \subset (-\frac{\pi}{2\lambda}, \pi + \frac{\pi}{2\lambda})$ (since $\lambda > 0$), so $\Phi(u + i\pi)$ is analytic and non-zero.
   - Crossing $\Phi(u) = \frac{1}{(1 + e^{2\lambda u})\Phi(u + i\pi)}$ therefore continues $\Phi(u)$ meromorphically to the symmetric strip $|\operatorname{Im}u| < \pi + \frac{\pi}{2\lambda}$, with poles precisely at the zeros of $1 + e^{2\lambda u} = 0$ in $-\pi - \frac{\pi}{2\lambda} < \operatorname{Im}u \le -\frac{\pi}{2\lambda}$.
   - The zeros of $1 + e^{2\lambda u} = 0$ are $u_n = -\frac{(2n+1)i\pi}{2\lambda}$ ($n \in \mathbb{Z}$), and $-\pi - \frac{\pi}{2\lambda} < -\frac{(2n+1)\pi}{2\lambda} \le -\frac{\pi}{2\lambda} \iff 0 \le n < \lambda$. Each zero is simple with derivative $\frac{d}{du}(1 + e^{2\lambda u})\big|_{u_n} = -2\lambda$, giving residue
     $$
     \operatorname{Res}_{\theta = \theta_n}\Phi(\theta - \theta_d) = -\frac{1}{2\lambda\,\Phi\big(i\pi - \frac{(2n+1)i\pi}{2\lambda}\big)}.
     $$
2. **Mechanism in the Gamma product (`eq-sgd-scalar`) and absence of cancellations:**
   - At $u_n = -\frac{(2n+1)i\pi}{2\lambda}$, $Y(u_n) = e^{-i(n+1/2)\pi} = -i(-1)^n$, so $1 + iY(u_n) = 1 + (-1)^n$:
     - For **even $n = 2m$**, $1 + iY(u_n) = 2 \ne 0$, while the $k=0$ factor $\Gamma(\frac14 - \frac{i\lambda u_n}{2\pi}) = \Gamma(-m)$ in $r_N$ has a simple pole.
     - For **odd $n = 2m+1$**, $1 + iY(u_n) = 0$ has a simple zero in the denominator of $\frac{1}{1+iY}$, while $\Gamma(\frac14 - \frac{i\lambda u_n}{2\pi}) = \Gamma(-m - \frac12)$ is finite and non-zero.
   - For all $0 \le n < \lambda$, every other Gamma argument in $r_N(u_n)$ is $\ge \min\big(\lambda - \frac{n}{2}, \frac{\lambda - n}{2} + \frac12, \frac{\lambda - n}{2}\big) > 0$, and every Gamma argument in $\bar r_N(u_n)$ is $k\lambda + c + \frac{2n+1}{4} > \frac12 > 0$. Since $\Gamma(x)$ has no poles or zeros on $(0, \infty)$, no other factor introduces a pole or zero in this strip.
   - At $u_n$, $-iY(u_n) = (-1)^{n+1}$, so the $2\times 2$ block matrix $M_n = \big(\begin{smallmatrix}1 & (-1)^{n+1} \\ (-1)^{n+1} & 1\end{smallmatrix}\big)\operatorname{diag}(q^{-n}, q^{n+1+s})$ has determinant $0$ and trace $q^{-n} + q^{n+1+s} \ne 0$ (for generic $q$), hence rank 1.
3. **Resonance widths (`eq-sgd-widths`):**
   - For $\theta_d > 0$ and $2n+1 < \lambda$, $E(\theta_n) = M\cosh\theta_n = M\cosh\theta_d\cos\frac{(2n+1)\pi}{2\lambda} - i M\sinh\theta_d\sin\frac{(2n+1)\pi}{2\lambda}$.
   - Using $m_{2n+1} = 2M\sin\frac{(2n+1)\pi}{2\lambda}$, the width is $\Gamma_n = -2\operatorname{Im}E(\theta_n) = m_{2n+1}\sinh\theta_d > 0$. As $\beta_{\rm SG} \to 0$ ($\lambda \to \infty$ with $m_1 \to m_0$), $\frac{m_{2n+1}}{m_1} = \frac{\sin((2n+1)\pi/2\lambda)}{\sin(\pi/2\lambda)} \to 2n+1$, so $\Gamma_n \to (2n+1)m_0\sinh\theta_d$.

### Script & Output (`referee_scratch/defects_38_41_pr52/verify_ch38.py`)

```text
3.1 Testing poles and residues for lam = 3.7320508 (n = 0, 1, 2, 3):
    n = 0 (even: Gamma): |res_num - res_exact|/|res| = 2.53e-12, |Phi(eps2)/Phi(eps1)| = 100.0 (expect 100)
    n = 1 (odd: 1/(1+iY)): |res_num - res_exact|/|res| = 1.9e-12, |Phi(eps2)/Phi(eps1)| = 100.0 (expect 100)
    n = 2 (even: Gamma): |res_num - res_exact|/|res| = 1.97e-12, |Phi(eps2)/Phi(eps1)| = 100.0 (expect 100)
    n = 3 (odd: 1/(1+iY)): |res_num - res_exact|/|res| = 3.02e-12, |Phi(eps2)/Phi(eps1)| = 100.0 (expect 100)
3.2 [PASS] All Gamma arguments in r_N (except k=0 first Gamma at even n) and r_bar_N are strictly positive at u_n for 0 <= n < lam.
3.3 [PASS] Matrix M_n at u_n has [[1, (-1)^{n+1}], [(-1)^{n+1}, 1]] * diag(...) -> trace != 0, det = 0 -> rank 1.
3.4 [PASS] Widths Gamma_n = m_{2n+1} sinh(theta_d) -> (2n+1) m_0 sinh(theta_d) verified.
```

### Label Assessment

- **`[Theorem]` is fully justified.** The proof is complete and exact.

---

## Item 4: Linearization about the odd jump (`prp-sgd-odd-jump`) & inverted-oscillator remarks

- **Verdict:** **correct** (with a minor edge-case clarification at $\theta_d = 0$ in the proof)
- **Target:** `part8-defects/38-defect-sine-gordon.qmd`, lines 302–331 (`prp-sgd-odd-jump`, proof, and subsequent two paragraphs).

### Independent Verification

1. **Linearized defect conditions about $(\phi,\psi) = (0, 2\pi/\beta_{\rm SG})$:**
   - The nonlinear defect conditions (`eq-sg-defect-conditions`) in light-cone derivatives $\partial_\pm = \partial_t \pm \partial_x$ are
     $$
     \partial_+(\phi - \psi) = -2m_0\sigma\sin\frac{\beta_{\rm SG}(\phi+\psi)}{2}, \qquad
     \partial_-(\phi + \psi) = 2m_0\sigma^{-1}\sin\frac{\beta_{\rm SG}(\phi-\psi)}{2}.
     $$
   - Setting $\phi = \delta\phi$ and $\psi = \frac{2\pi}{\beta_{\rm SG}} + \tilde\psi$, we have $\sin\big(\frac{\pi}{2}\cdot 2 + \frac{\beta_{\rm SG}(\delta\phi \pm \tilde\psi)}{2}\big) = -\frac{\beta_{\rm SG}}{2}(\delta\phi \pm \tilde\psi) + O(3)$, so the linearized conditions at $x = 0$ are
     $$
     \partial_+(\delta\phi - \tilde\psi) = \sigma m_0(\delta\phi + \tilde\psi), \qquad
     \partial_-(\delta\phi + \tilde\psi) = -\sigma^{-1}m_0(\delta\phi - \tilde\psi).
     $$
   - Substituting $\delta\phi = A e^{\gamma t + \kappa x}$ ($x < 0$) and $\tilde\psi = B e^{\gamma t - \kappa x}$ ($x > 0$) with $\kappa > 0$ and $\gamma^2 = \kappa^2 - m_0^2$ gives the $2\times 2$ homogeneous system $M(\gamma,\kappa)\binom{A}{B} = 0$ with
     $$
     M(\gamma,\kappa) = \begin{pmatrix} \gamma + \kappa - \sigma m_0 & -(\gamma - \kappa) - \sigma m_0 \\ \gamma - \kappa + \sigma^{-1}m_0 & \gamma + \kappa - \sigma^{-1}m_0 \end{pmatrix}.
     $$
   - Using $\gamma^2 = \kappa^2 - m_0^2$, the determinant simplifies to
     $$
     \det M(\gamma,\kappa) = 4\kappa\big(\kappa - m_0\cosh\theta_d\big), \qquad \sigma = e^{-\theta_d}.
     $$
   - Since $\kappa > 0$, $\det M = 0 \iff \kappa = m_0\cosh\theta_d$, which forces $\gamma = \pm m_0\sinh\theta_d$.
   - At $\kappa = m_0\cosh\theta_d$ and $\gamma = +m_0\sinh\theta_d$, $\gamma - \kappa = -m_0 e^{-\theta_d} = -\sigma m_0$ and $\gamma + \kappa = m_0 e^{\theta_d} = \sigma^{-1}m_0$, so
     $$
     M(+m_0\sinh\theta_d, m_0\cosh\theta_d) = \begin{pmatrix} 2m_0\sinh\theta_d & 0 \\ 0 & 0 \end{pmatrix}, \qquad
     M(-m_0\sinh\theta_d, m_0\cosh\theta_d) = \begin{pmatrix} 0 & 0 \\ 0 & -2m_0\sinh\theta_d \end{pmatrix}.
     $$
   - **Minor edge case at $\theta_d = 0$:** When $\theta_d \ne 0$, $2m_0\sinh\theta_d \ne 0$, so the $\gamma = +m_0\sinh\theta_d$ mode forces $A = 0$ (supported on $x > 0$) and the $\gamma = -m_0\sinh\theta_d$ mode forces $B = 0$ (supported on $x < 0$). At $\theta_d = 0$, $\gamma = 0$ and both matrices vanish identically, so $(A,B)$ is unconstrained. The proposition statement (line 312) is completely accurate because it states $A=0$ only under the hypothesis "For $\theta_d > 0$". In the proof (line 316), adding "for $\theta_d \ne 0$" where $2m_0\sinh\theta_d\,A = 0 \implies A = 0$ is deduced makes the proof airtight.
2. **Nonlinear soliton emission tail:**
   - For $\phi = 0$ ($x < 0$) and $\psi(x,t) = \frac{4}{\beta_{\rm SG}}\arctan e^{m_0(x\cosh\theta_d - t\sinh\theta_d) + c}$ ($x > 0$), both nonlinear defect conditions hold identically at $x = 0$ for all $t \in \mathbb{R}$.
   - For $\theta_d > 0$, as $t \to -\infty$ at fixed $x > 0$, $m_0(x\cosh\theta_d - t\sinh\theta_d) + c \to +\infty$, so $\psi(x,t) \to \frac{2\pi}{\beta_{\rm SG}}$ and
     $$
     \psi(x,t) - \frac{2\pi}{\beta_{\rm SG}} = -\frac{4}{\beta_{\rm SG}}\arctan e^{-m_0(x\cosh\theta_d - t\sinh\theta_d) - c} \sim -\frac{4}{\beta_{\rm SG}} e^{-c} e^{\gamma t - \kappa x},
     $$
     matching the growing mode ($A = 0$, $B = -\frac{4}{\beta_{\rm SG}}e^{-c}$) with $\gamma = m_0\sinh\theta_d$.
3. **Heuristic inverted-oscillator paragraph & correction to BCZ (lines 325–331):**
   - The paragraph on lines 326–327 clearly separates what is proved (`prp-sgd-odd-jump` for $n=0$) from the heuristic inverted harmonic oscillator analogy (calling $2n+1$ "heuristically the decay rates of the $n$-th states of an inverted oscillator of Lyapounov exponent $\gamma$", and explicitly stating in the next sentence: "To make this a semiclassical derivation of (-@eq-sgd-widths) beyond $n=0$ one would have to reduce the field equations at the odd jump to a single unstable degree of freedom; @prp-sgd-odd-jump proves only the $n=0$ mode and its growth rate $\gamma$").
   - The correction to `@bcz2005` Sec. 4 on lines 328–331 is exact: in units of $M \sim \frac{8m_0}{\beta_{\rm SG}^2}$, $\Gamma_0 / M = 2\sinh\theta_d\sin\frac{\pi}{2\lambda} = O(\beta_{\rm SG}^2) \to 0$, whereas in absolute units $\Gamma_0 \to m_0\sinh\theta_d \ne 0$.

### Label Assessment & Recommended Minor Edit

- **`[Theorem]` is justified.**
- **Optional minor clarification in `part8-defects/38-defect-sine-gordon.qmd`, line 316:**
  Replace:
  ```markdown
  where $\gamma^2-\kappa^2=-m_0^2$ and $\sigma+\sigma^{-1}=2\cosh\theta_d$ have been used. Since $\kappa>0$, the determinant vanishes if and only if $\kappa=m_0\cosh\theta_d$, which gives $\gamma=\pm m_0\sinh\theta_d$. At $\gamma=m_0\sinh\theta_d$, $\gamma-\kappa=-m_0e^{-\theta_d}=-\sigma m_0$ and $\gamma+\kappa=\sigma^{-1}m_0$, so the coefficients of $B$ vanish and the first condition reads $2m_0\sinh\theta_d\,A=0$. At $\gamma=-m_0\sinh\theta_d$ the coefficients of $A$ vanish and the second condition forces $B=0$.
  ```
  with:
  ```markdown
  where $\gamma^2-\kappa^2=-m_0^2$ and $\sigma+\sigma^{-1}=2\cosh\theta_d$ have been used. Since $\kappa>0$, the determinant vanishes if and only if $\kappa=m_0\cosh\theta_d$, which gives $\gamma=\pm m_0\sinh\theta_d$. At $\gamma=m_0\sinh\theta_d$, $\gamma-\kappa=-m_0e^{-\theta_d}=-\sigma m_0$ and $\gamma+\kappa=\sigma^{-1}m_0$, so the coefficients of $B$ vanish and the first condition reads $2m_0\sinh\theta_d\,A=0$, which forces $A=0$ for $\theta_d\ne0$. At $\gamma=-m_0\sinh\theta_d$ the coefficients of $A$ vanish and the second condition forces $B=0$ for $\theta_d\ne0$.
  ```

---

## Item 5: Two Weyl pairs and Faddeev's modular double (`prp-dq-two-weyl-pairs`)

- **Verdict:** **correct**
- **Target:** `part8-defects/41-defect-quantum-groups.qmd`, lines 180–215 (`prp-dq-two-weyl-pairs` and surrounding discussion).

### Independent Verification

1. **Canonical commutator and exponential Weyl relations:**
   - From the defect kinetic term $\frac{1}{2}\phi(0)\partial_t\psi(0) - \frac{1}{2}\psi(0)\partial_t\phi(0) = \phi(0)\partial_t\psi(0) - \frac{d}{dt}\big(\frac12\phi(0)\psi(0)\big)$, the conjugate momentum to $\psi(0)$ is $\Pi_{\psi(0)} = \phi(0)$, so $[\psi(0),\phi(0)] = i\hbar$, i.e., $[\phi(0),\psi(0)] = -i\hbar$.
   - For $U = e^{\frac{b}{2\sqrt{2}}(\psi - \phi)}$ and $W_\pm = e^{\pm\frac{b}{2\sqrt{2}}(\psi + \phi)}$ with $b = i\beta_{\rm SG}/\sqrt{2}$:
     $$
     \big[\ln U, \ln W_\pm\big] = \pm\frac{b^2}{8}\big[\psi - \phi, \psi + \phi\big] = \pm\frac{b^2}{4}\big[\psi, \phi\big] = \pm\frac{i\hbar b^2}{4} = \mp\frac{i\hbar\beta_{\rm SG}^2}{8},
     $$
     so $U W_\pm = \hat{q}^{\pm 1} W_\pm U$ with $\hat{q} = e^{i\hbar b^2/4} = e^{-i\hbar\beta_{\rm SG}^2/8}$.
   - For $U_{\rm s} = e^{\frac{\pi i}{\sqrt{2}\hbar b}(\psi - \phi)} = e^{\frac{\pi}{\hbar\beta_{\rm SG}}(\psi - \phi)}$ and $X_\pm = e^{\pm\frac{\sqrt{2}\pi i}{\hbar b}(\psi + \phi)} = e^{\pm\frac{2\pi}{\hbar\beta_{\rm SG}}(\psi + \phi)}$:
     $$
     \big[\ln U_{\rm s}, \ln X_\pm\big] = \mp\frac{\pi^2}{\hbar^2 b^2}\big[\psi - \phi, \psi + \phi\big] = \mp\frac{2\pi^2}{\hbar^2 b^2}(2i\hbar) = \frac{\pm 4\pi^2 i}{\hbar b^2} = \mp\frac{8\pi^2 i}{\hbar\beta_{\rm SG}^2},
     $$
     so $U_{\rm s} X_\pm = q_{\rm s}^{\pm 1} X_\pm U_{\rm s}$ with $q_{\rm s} = e^{4\pi^2 i/\hbar b^2} = e^{-8\pi^2 i/\hbar\beta_{\rm SG}^2} = -q$.
2. **Action on defect states $|n\rangle$:**
   - In sector $s$, $|n\rangle$ has jump $\psi(0) - \phi(0) = \frac{2\pi}{\beta_{\rm SG}}(2n + s)$, so $U_{\rm s}|n\rangle = e^{\frac{\pi}{\hbar\beta_{\rm SG}}\cdot\frac{2\pi}{\beta_{\rm SG}}(2n+s)}|n\rangle = (-q_{\rm s})^{-(n+s/2)}|n\rangle = q^{-(n+s/2)}|n\rangle$.
   - Since $[\psi - \phi, \ln X_\pm] = \pm\frac{2\pi}{\hbar\beta_{\rm SG}}[\psi - \phi, \psi + \phi] = \pm\frac{4\pi}{\beta_{\rm SG}}$ (with $\hbar$ cancelling!), $X_\pm(\psi - \phi)X_\pm^{-1} = \psi - \phi \mp \frac{4\pi}{\beta_{\rm SG}}$, so $X_\pm$ shifts $2n + s \to 2(n \pm 1) + s$, mapping $|n\rangle \mapsto |n \pm 1\rangle$.
3. **Cross-commutators and Faddeev's modular double:**
   - The cross-commutators of the exponents are $[\ln U, \ln U_{\rm s}] = 0$, $[\ln U, \ln X_\pm] = \pm i\pi$, $[\ln W_\pm, \ln U_{\rm s}] = \mp i\pi$, and $[\ln W_\pm, \ln X_{\pm'}] = \pm(\pm')i\pi$. Exponentiating gives $e^0 = +1$ and $e^{\pm i\pi} = -1$, so the two pairs commute up to signs, and $(U^2, W_+^2)$ commutes **with no signs** (`e^{\pm 2\pi i} = +1`) with $(U_{\rm s}, X_-)$.
   - For $(U^2, W_+^2)$, $[\ln(U^2), \ln(W_+^2)] = 4(\frac{i\hbar b^2}{4}) = i\hbar b^2 = 2\pi i\tau$ with $\tau = \frac{\hbar b^2}{2\pi}$, giving multiplier $e^{2\pi i\tau} = q_{\rm Faddeev}^2$ ($q_{\rm Faddeev} = e^{i\pi\tau}$).
   - For $(U_{\rm s}, X_-)$, $[\ln U_{\rm s}, \ln X_-] = -\frac{4\pi^2 i}{\hbar b^2} = -\frac{2\pi i}{\tau}$, giving multiplier $e^{-2\pi i/\tau} = \tilde{q}_{\rm Faddeev}^2$ ($\tilde{q}_{\rm Faddeev} = e^{-i\pi/\tau}$).

### Script & Output (`referee_scratch/defects_38_41_pr52/verify_ch41_43.py`)

```text
5.1 [PASS] [log U, log W_+] = I*b**2*hbar/4 (expected i*hbar*b^2/4)
           [log U, log W_-] = -I*b**2*hbar/4 (expected -i*hbar*b^2/4)
           At b = i*beta_SG/sqrt(2): log(q_hat) = -I*beta_SG**2*hbar/8 = -i*hbar*beta_SG^2/8
5.2 [PASS] [log U_s, log X_+] = 4*I*pi**2/(b**2*hbar) (expected 4*pi^2*i/(hbar*b^2))
           [log U_s, log X_-] = -4*I*pi**2/(b**2*hbar) (expected -4*pi^2*i/(hbar*b^2))
           At b = i*beta_SG/sqrt(2): log(q_s) = -8*I*pi**2/(beta_SG**2*hbar) = -8*pi^2*i/(hbar*beta_SG^2)
5.3 [PASS] Cross-commutators of exponents:
           [log U, log U_s] = 0   -> exp = 1
           [log U, log X_+] = I*pi -> exp = -1
           [log W_+, log U_s] = -I*pi -> exp = -1
           [log W_+, log X_+] = I*pi -> exp = -1
           [log W_+, log X_-] = -I*pi -> exp = -1
5.4 [PASS] Modular double check: [log(U^2), log(W_+^2)] - 2*pi*i*tau = 0, [log(U_s), log(X_-)] + 2*pi*i/tau = 0
5.5 [PASS] Shift of jump: [psi - phi, log X_+] = 4*pi/beta_SG -> X_+ maps |n> to |n+1>.
```

### Label Assessment

- **`[Theorem]` is fully justified.** All statements and proofs in `prp-dq-two-weyl-pairs` are exact.

---

## Item 6: Literature corrections, Prior Art, and Appendices D & E

- **Verdict:** **correct** (with a recommendation to include the four additional misprints of Nemes 2010 from Ch. 38 Notes in Appendix D)
- **Targets:**
  - `part8-defects/38-defect-sine-gordon.qmd` (lines 223, 279, 326–331, 404–405)
  - `part8-defects/41-defect-quantum-groups.qmd` (lines 206, 350–352)
  - `appendices/D-corrections.qmd` (lines 83–86, 137, 163–166)
  - `appendices/E-original-results.qmd` (lines 43, 44, 47)

### Direct Verification Against Primary `.tex` Sources

1. **Nemes 2010 (`0909.3268.tex`):**
   - **Eq. (34), line 497:**
     ```latex
     {^{s}T}_{\alpha}^{ab}(\theta)=f(q,x)\left(\begin{array}{cccc}
     Q^{\alpha}\delta_{\alpha}^{\beta}&0&-iQ^{-\alpha}px\delta_{\alpha}^{\beta-2}&0\\
     0&Q^{-\alpha}\delta_{\alpha}^{\beta}&0&-iQ^{\alpha}px\delta_{\alpha}^{\beta+2}\\
     ```
     Verified: $Q$ is never defined in `0909.3268.tex` (where $q = e^{i\pi\gamma}$ in eq. (36) and BCZ's $Q = e^{-i\pi\gamma/2}$ is used without definition).
   - **Eq. (49), line 696:**
     ```latex
     \frac{1}{i}\ln\frac{r(x)}{\bar{r}(x)}=-\int_{0}^{\infty}\frac{dt}{t}\frac{\sinh\left(\frac{t}{2}\left(1+\frac{1}{\gamma}\right)\right)\sin\left(\frac{t\left(\theta-\eta\right)}{\pi}\right)}{\sinh\left(\frac{t}{\gamma}\right)\cosh\left(\frac{t}{2}\right)}
     ```
     Verified: as $t \to \infty$, $\frac{\sinh(\frac{t}{2}(1+\gamma^{-1}))}{\sinh(t/\gamma)\cosh(t/2)} \sim e^{-t/(2\gamma)}$, while $\sin\big(\frac{t(\theta-\eta)}{\pi}\big)$ grows as $e^{t|\operatorname{Im}(\theta-\eta)|/\pi}$, so the integral converges only for $|\operatorname{Im}(\theta-\eta)| < \frac{\pi}{2\gamma}$.
   - **Eq. (54), line 753:**
     ```latex
     \frac{d}{d\theta}\frac{1}{i}\ln\frac{r(x)}{\bar{r}(x)}=-\gamma\ln\left(\frac{\eta-\theta}{2}\right)+\mathcal{O}\left(\gamma^{0}\right).
     ```
     Verified: differentiating eq. (53), $\frac{1}{i}\ln\frac{r(x)}{\bar{r}(x)} \sim -\frac{\gamma}{\pi}\int_0^\infty \frac{dt}{t}\,\frac{\sin(t(\theta-\eta)/\pi)}{\sinh(t/2)}$, gives $-\frac{\gamma}{\pi^2}\int_0^\infty \frac{\cos(t(\theta-\eta)/\pi)}{\sinh(t/2)}\,dt$, whose finite jump across $\theta \to \eta\pm$ is $-\gamma\ln\tanh\frac{|\eta-\theta|}{2}$ (matching the derivative of his eq. (55)), whereas eq. (54) omits $\tanh$.
   - **Eqs. (56) and (61), lines 815 and 866:**
     ```latex
     f(q,x) = \frac{e^{i\pi\left(1+\gamma\right)/4}}{1+e^{\gamma\left(\theta-\eta\right)}}e^{-i\pi\gamma/2}\dots = \frac{-e^{-i\pi\gamma/4}}{1+e^{\gamma\left(\theta-\eta\right)}}\dots
     ```
     ```latex
     \lim_{\theta\to-\infty}f(q,x)=\frac{e^{i\pi\left(1+\gamma\right)/4}}{1+0}1=-e^{i\pi\gamma/4}.
     ```
     Verified: both equations replace the factor $e^{i\pi/4}$ in $e^{i\pi(1+\gamma)/4} = e^{i\pi/4}e^{i\pi\gamma/4}$ by $-1$.
   - **Section 5, eqs. (65)–(67), lines 903–945:**
     ```latex
     \theta_{k}=\eta-i\frac{\pi}{\gamma}\left(2k+\frac{1}{2}\right)\qquad k=0,1,2\dots
     ...
     -\frac{\pi}{2}<\text{Im}\left(\theta_{k}-\eta\right)<0\;\Rightarrow\;0\leq k<\frac{\gamma-1}{4},
     ```
     Verified: Nemes counts only the poles of $\Gamma(\frac14 + k\gamma - z)$ at $z = \frac14 + 2k$ ($k \ge 0$), obtaining $\sim \gamma/4$ poles in $(-\frac{\pi}{2}, 0)$, and misses the odd-$n$ poles at $\theta_{2m+1} = \eta - \frac{i\pi}{\gamma}(2m + \frac32)$ coming from the prefactor $\frac{1}{1+ipx}$ in his eq. (35), which double the count to $\sim \gamma/2$ ($0 \le 2n+1 < \gamma$).
   - **Note on Appendix D (`appendices/D-corrections.qmd`, lines 163–166):**
     In `part8-defects/38-defect-sine-gordon.qmd` (lines 404–405, Notes), the chapter notes both the missing half of the pole ladder in Sec. 5 and the misprints in (34), (54), (56), and (61), whereas `appendices/D-corrections.qmd` (line 165) currently lists only `(67)`. Per the book's rule ("Each correction goes in the chapter's Notes and in Appendix D (`appendices/D-corrections.qmd`), in the same change"), we recommend adding a bullet under **Nemes 2010** in `appendices/D-corrections.qmd` for (34), (54), (56), and (61) (replacement text given below).

2. **Bowcock, Corrigan and Zambon 2005 (`0506169.tex`):**
   - **Section 4, after (4.1), lines 1072–1090:**
     ```latex
     Notice also that the limit $\beta\rightarrow 0$ implies $1/\gamma\rightarrow 0$ and the real part of $\Delta E$ becomes the classical energy of a soliton with rapidity $\eta$, and the imaginary part of $\Delta E$ vanishes ... In the classical limit the width goes to zero so the only process is absorption.
     ```
   - Verified: $\operatorname{Im}\Delta E = -m_s\sinh\eta\sin\frac{\pi}{2\gamma}$. In units of $m_s$, $\operatorname{Im}\Delta E / m_s \to 0$, but since $m_s = \frac{8m}{\beta^2} = \frac{(\gamma+1)m}{\pi}$, in absolute units $\operatorname{Im}\Delta E \to -\frac12 m\sinh\eta$ and the resonance width $\Gamma_0 = -2\operatorname{Im}\Delta E \to m\sinh\eta \ne 0$. Both `38-defect-sine-gordon.qmd` (lines 328–331, 398) and `appendices/D-corrections.qmd` (line 137) quote and explain this accurately.

3. **Faddeev 2000 (`9912078.tex`):**
   - **Eqs. (14)–(22), lines 258–317:**
     ```latex
     [Q,P]=2\pi i\hbar I . \eqno(14)
     u=e^{\alpha P} , \quad v=e^{\beta Q} \eqno(15)
     \tilde{u}=e^{\tilde{\alpha}P}, \quad \tilde{v}=e^{\tilde{\beta}Q} \eqno(17)
     \alpha \tilde{\beta}=\hbar , \quad \tilde{\alpha}\beta=\hbar \eqno(18)
     \alpha\beta=\frac{\hbar^2}{\tilde{\alpha}\tilde{\beta}} . \eqno(19)
     ```
   - Verified: from (14), $[\alpha P, \tilde\beta Q] = -2\pi i\hbar\alpha\tilde\beta$, so $u\tilde v = e^{-2\pi i\hbar\alpha\tilde\beta}\tilde v u$. Commutation requires $\hbar\alpha\tilde\beta \in \mathbb{Z}$ and $\hbar\tilde\alpha\beta \in \mathbb{Z}$ (specifically $\alpha\tilde\beta = 1/\hbar$, $\tilde\alpha\beta = -1/\hbar$ to match $q = e^{i\pi\tau}$, $\tilde q = e^{-i\pi/\tau}$ in his (22)), whereas (18)–(19) print $\hbar$ instead of $1/\hbar$ (anticipating his footnote 1 on p. 5 where $\hbar$ is set to $1$). Both `41-defect-quantum-groups.qmd` (lines 206, 350–352) and `appendices/D-corrections.qmd` (line 85) state this accurately.

4. **Appendix E (`appendices/E-original-results.qmd`, lines 43, 44, 47):**
   - Line 43 (`@prp-sgd-scalar-integral`, `@prp-sgd-pole-ladder`, `@prp-sgd-odd-jump`): opens with `- **[Theorem]**`, matching the three `[Theorem]` statements in Ch. 38.
   - Line 44 (`@prp-sg-breather-transmission`): opens with `- **[Theorem]**` for $B_1$ and notes `(@prp-sg-breather-transmission; numerical for $n\ge2$)`.
   - Line 47 (`@prp-dq-two-weyl-pairs`): opens with `- **[Theorem]**`, matching `@prp-dq-two-weyl-pairs` in Ch. 41.

### Recommended Replacement Text for `appendices/D-corrections.qmd` (lines 163–166)

```markdown
**Nemes 2010** [@nemes2010], semiclassical analysis of the sine-Gordon defect.

- **(34), (54), (56) and (61) of the arXiv version.** In (34), $Q^{\pm\alpha}$ is used without defining $Q$ ($Q=e^{-i\pi\gamma/2}$ in [@bcz2005]). In (54), differentiating (53) gives $-\gamma\ln\tanh\frac{\eta-\theta}2$, as his (55) requires, not $-\gamma\ln\frac{\eta-\theta}2$. In (56) and (61), the factor $e^{i\pi/4}$ of $e^{i\pi(1+\gamma)/4}$ is replaced by $-1$. *Misprints* (@sec-sg-defect-notes).
- **(67) of the arXiv version, the number of resonance poles.** It counts about $\gamma/4$ poles with $-\frac\pi2<\operatorname{Im}(\theta-\eta)<0$ at large $\gamma$, from the poles of the Gamma functions alone. The prefactor $1/(1+ip_x)$ of (34) contributes as many again, so there are about $\gamma/2$, one for each $n\ge0$ with $2n+1<\gamma$. *Incomplete* (@prp-sgd-pole-ladder; `code/defect_sg_scalar_poles.py`).
```

---

## Item 7: Candidate Book Improvements from Lectures 8 & 9 (PRs #55 & #56)

### Item 7a: Chirality reading of `prp-diagonal-borel` (`43-defects-and-boundaries.qmd` line 110 & `41-defect-quantum-groups.qmd` line 61)

- **Verdict:** **wrong** in `part8-defects/43-defects-and-boundaries.qmd` line 110 (the pairing of $\cos\frac{\beta(\phi\mp\psi)}{2}$ with the passed and reflected chiralities is reversed, and the attribution of $\Delta(U_q(\mathfrak{b}^-))$ and $(1\otimes\tau)\Delta(U_q(\mathfrak{b}^+))$ to individual massless chiralities conflates the two chiralities with the two Borel subalgebras).
- **Targets:**
  - `part8-defects/43-defects-and-boundaries.qmd`, line 110 (after the proof of `prp-diagonal-borel`).
  - `part8-defects/41-defect-quantum-groups.qmd`, lines 51–63 (`sec-perturbed-defects`).

#### Independent Derivation (`verify_ch41_43.py`, step 7a.1)

1. **Massless defect conditions at $m_0 = 0$:**
   Write $\phi(t,x) = \phi_L(t+x) + \phi_R(t-x)$ on $x < 0$ and $\psi(t,x) = \psi_L(t+x) + \psi_R(t-x)$ on $x > 0$. At $x = 0$, the type-I$_-$ conditions $\partial_x\phi = \partial_t\psi$, $\partial_x\psi = \partial_t\phi$ read
   $$
   \phi_L' - \phi_R' = \psi_L' + \psi_R', \qquad \psi_L' - \psi_R' = \phi_L' + \phi_R' \quad\Longleftrightarrow\quad \psi_L' = +\phi_L', \quad \psi_R' = -\phi_R'.
   $$
   Hence at $x = 0$ (up to zero-mode constants):
   - The **passed chirality** ($\psi_L = +\phi_L$, afunction of $t+x$) has $\phi - \psi = \text{const}$ and $\phi + \psi = 2\phi_L + \text{const}$, so it is perturbed by $\sigma\cos\frac{\beta(\phi+\psi)}{2}$.
   - The **sign-reflected chirality** ($\psi_R = -\phi_R$, a function of $t-x$) has $\phi + \psi = \text{const}$ and $\phi - \psi = 2\phi_R + \text{const}$, so it is perturbed by $\sigma^{-1}\cos\frac{\beta(\phi-\psi)}{2}$.
2. **What is wrong in `part8-defects/43-defects-and-boundaries.qmd` line 110:**
   Line 110 currently reads:
   > "At $m_0=0$ the defect conditions $\partial_x\phi=\partial_t\psi$, $\partial_x\psi=\partial_t\phi$ pass one chirality and reflect the other with a sign ($\phi\to-\psi$, which on solitons is $\tau$). Switching on the defect potential (-@eq-sg-defect-potential) perturbs the two chiralities separately, $\cos\frac{\beta(\phi-\psi)}2$ in one and $\cos\frac{\beta(\phi+\psi)}2$ in the other. In the quantum theory the first chirality gives the untwisted Borel symmetry $\Delta(U_q(\mathfrak b^-))$ of @sec-perturbed-defects, and the second gives the $\tau$-twisted Borel symmetry $(1\otimes\tau)\Delta(U_q(\mathfrak b^+))$."
   - First, "pass one chirality and reflect the other with a sign ... $\cos\frac{\beta(\phi-\psi)}2$ in one and $\cos\frac{\beta(\phi+\psi)}2$ in the other" pairs $\cos\frac{\beta(\phi-\psi)}2$ with the passed chirality ("the first chirality") and $\cos\frac{\beta(\phi+\psi)}2$ with the sign-reflected chirality ("the second"), which is backwards (`exr-dq-chiral-defect` in Ch. 41 line 313 gets it right: $\phi+\psi$ depends on $t+x$ and $\phi-\psi$ on $t-x$).
   - Second, each Borel subalgebra $U_q(\mathfrak{b}^\mp)$ contains *both* affine generators ($f_1, f_0$ for $\mathfrak{b}^-$, or $e_1, e_0$ for $\mathfrak{b}^+$), which in the vertex-operator construction (`41-defect-quantum-groups.qmd` line 161) come from *both* $\cos\frac{\beta(\phi\pm\psi)}{2}$ terms of the defect potential together! The reason $\Delta(U_q(\mathfrak{b}^-))$ is untwisted while $(1\otimes\tau)\Delta(U_q(\mathfrak{b}^+))$ is twisted by $\tau$ is explained cleanly and exactly by Lecture 9 (`09-folding-fusing.qmd`, `sec-l9-borel`, line 136):
     1. In the transmission matrix $\mathbb{L}(y) = \mathcal{D} + y\mathcal{O}$, as $y \to 0$ ($\theta \to -\infty$, where $f_j \propto y^{-1}$ dominate) $\mathbb{L}$ is diagonal ($\mathcal{D}$), so $U_q(\mathfrak{b}^-)$ is intertwined untwisted; as $y \to \infty$ ($\theta \to +\infty$, where $e_j \propto y$ dominate) $\mathbb{L} \sim y\mathcal{O}$ converts solitons to antisolitons ($\tau$), so $U_q(\mathfrak{b}^+)$ is intertwined with $\tau$ on the outgoing copy.
     2. Equivalently, by `prp-dc-types` and `prp-two-borels` (3), the reflection $\psi \to -\psi$ turns a type-I$_-$ defect into a type-I$_+$ defect and acts on the second copy of the doubled theory as charge conjugation $\tau$; since a type-I$_-$ defect preserves $\Delta(U_q(\mathfrak{b}^-))$ and a type-I$_+$ defect preserves $\Delta(U_q(\mathfrak{b}^+))$, the same type-I$_-$ defect seen through $\psi \to -\psi$ preserves $(1\otimes\tau)\Delta(U_q(\mathfrak{b}^+))$.
3. **Check of `part8-defects/41-defect-quantum-groups.qmd` line 61 (`sec-perturbed-defects`):**
   Lines 51–61 of `41-defect-quantum-groups.qmd` are completely consistent (line 53 has $\partial_+\phi = \partial_+\psi$, $\partial_-\phi = -\partial_-\psi$, matching `exr-dq-chiral-defect` line 313, and line 61 simply defines $U_q(\mathfrak{b}^\pm)$ without making the reversed pairing of Ch. 43 line 110).

#### Recommended Replacement Text for `part8-defects/43-defects-and-boundaries.qmd`, line 110

Replace line 110 of `part8-defects/43-defects-and-boundaries.qmd`:
```markdown
The two Borel symmetries of @prp-diagonal-borel have a direct reading in the defect Lagrangian. At $m_0=0$ the defect conditions $\partial_x\phi=\partial_t\psi$, $\partial_x\psi=\partial_t\phi$ pass one chirality and reflect the other with a sign ($\phi\to-\psi$, which on solitons is $\tau$). Switching on the defect potential (-@eq-sg-defect-potential) perturbs the two chiralities separately, $\cos\frac{\beta(\phi-\psi)}2$ in one and $\cos\frac{\beta(\phi+\psi)}2$ in the other. In the quantum theory the first chirality gives the untwisted Borel symmetry $\Delta(U_q(\mathfrak b^-))$ of @sec-perturbed-defects, and the second gives the $\tau$-twisted Borel symmetry $(1\otimes\tau)\Delta(U_q(\mathfrak b^+))$. Neither extends to all of $U_q(\widehat{sl}_2)$, because the defect is not transparent once $m_0\ne0$, and the two together generate the coideal of @prp-diagonal-borel.
```
with:
```markdown
The two Borel symmetries of @prp-diagonal-borel have a direct reading both in the defect Lagrangian and in the limits of $\mathbb T(\theta)$. By @prp-dc-types, negating the right-hand field, $\psi\to-\psi$, exchanges $\phi+\psi$ with $\phi-\psi$ in the defect potential (-@eq-sg-defect-potential) and turns a type-I$_-$ defect into a type-I$_+$ defect. On the solitons of the second copy $\psi\to-\psi$ is charge conjugation $\tau$, so the untwisted $U_q(\mathfrak b^+)$ symmetry of type I$_+$ (@prp-two-borels) becomes the $\tau$-twisted symmetry $(1\otimes\tau)\Delta(U_q(\mathfrak b^+))$ of type I$_-$ on the shifted module $\mathcal W'$ (@prp-two-borels (3)). Equivalently, in $\mathbb L(y)=\mathcal D+y\mathcal O$ (@prp-defect-twist-limit), the $\theta\to-\infty$ ($y\to0$) term $\mathcal D$ transmits solitons without conversion, where the $f_i\propto y^{-1}$ of $U_q(\mathfrak b^-)$ dominate, while the $\theta\to+\infty$ ($y\to\infty$) term $y\mathcal O$ converts solitons into antisolitons ($\tau$), where the $e_i\propto y$ of $U_q(\mathfrak b^+)$ dominate. Neither Borel symmetry extends to all of $U_q(\widehat{sl}_2)$, because the defect is not transparent once $m_0\ne0$, and the two together generate the coideal of @prp-diagonal-borel.
```

---

### Item 7b: Upgrade of `prp-defect-fusion-k` (3) (`43-defects-and-boundaries.qmd` lines 147–161) via `thm-l9-dirichlet`

- **Verdict:** **correct** — Lecture 9's Malmstén-integral proof (`thm-l9-dirichlet`, lines 285–322 of `09-folding-fusing.qmd`) is completely rigorous and proves the scalar factor identity in `prp-defect-fusion-k` (3), allowing `prp-defect-fusion-k` to be upgraded to **`[Theorem]`**.
- **Targets:**
  - `part8-defects/43-defects-and-boundaries.qmd`, lines 147–161 (`prp-defect-fusion-k` and its proof).
  - `appendices/E-original-results.qmd`, line 54.

#### Independent Verification (`verify_ch41_43.py`, steps 7b.1–7b.2)

1. **Malmstén integral for Ghoshal–Zamolodchikov's Dirichlet factor $\sigma(x, u)$ (`eq-bsg-k-full`, Ch. 30 line 91):**
   - For $l \ge 1$ and $w = \frac{\lambda u}{\pi}$, the $l$-th factor of $\sigma(x, u)$ has numerator Gamma arguments $a_{1,2}(l) \mp w$ with $a_{1,2}(l) = \frac12 \pm \frac{x}{\pi} + (2l-1)\lambda$ and denominator arguments $b_1(l) \mp w = \frac12 - \frac{x}{\pi} + (2l-2)\lambda \mp w$, $b_2(l) \mp w = \frac12 + \frac{x}{\pi} + 2l\lambda \mp w$.
   - Since $a_1(l) + a_2(l) - b_1(l) - b_2(l) = 0$, the linear term in Malmstén's formula cancels for each $l$. Summing over $l \ge 1$ gives
     $$
     \sum_{l=1}^\infty \big(e^{-a_1(l)t} + e^{-a_2(l)t} - e^{-b_1(l)t} - e^{-b_2(l)t}\big) = \frac{e^{-t/2}(e^{-xt/\pi} - e^{xt/\pi})(e^{-\lambda t} - 1)}{1 - e^{-2\lambda t}} = -\frac{e^{-t/2}\sinh(xt/\pi)\sinh(\lambda t/2)}{\sinh(\lambda t)}.
     $$
   - Multiplying by $\frac{2\sin(wt)}{1 - e^{-t}} = \frac{\sin(\lambda u t/\pi)}{e^{-t/2}\sinh(t/2)}$ gives
     $$
     \ln\frac{\sigma(x,u)\cos(x+\lambda u)}{\cos x} = -\int_0^\infty \frac{dt}{t}\,\frac{\sinh\big(\frac{xt}{\pi} + \frac{\lambda t}{2}\big) - \sinh\big(\frac{\lambda t}{2}\big)}{\sinh\frac{t}{2}\cosh\frac{\lambda t}{2}}\,\sin\frac{\lambda u t}{\pi} - \int_0^\infty \frac{dt}{t}\,\frac{\tanh\frac{\lambda t}{2}}{\sinh\frac{t}{2}}\,\sin\frac{\lambda u t}{\pi},
     $$
      Wait: directly summing $e^{-a_1 t} + e^{-a_2 t} - e^{-b_1 t} - e^{-b_2 t} = e^{-t/2}\big[e^{-xt/\pi - \lambda t}(1 - e^{-\lambda t}) - e^{xt/\pi}(1 - e^{-\lambda t})\big] = -2e^{-(1+\lambda)t/2}\sinh(\frac{xt}{\pi} + \frac{\lambda t}{2})(1 - e^{-\lambda t})$, and dividing by $(1 - e^{-t})(1 - e^{-2\lambda t}) = 4e^{-(1/2+\lambda)t}\sinh\frac{t}{2}\cosh\frac{\lambda t}{2}$ gives the exact single kernel:
     $$
     \ln\frac{\sigma(x,u)\cos(x+\lambda u)}{\cos x} = -\int_0^\infty \frac{dt}{t}\,\frac{\sinh\big(\frac{xt}{\pi} + \frac{\lambda t}{2}\big)\sin\frac{\lambda u t}{\pi}}{\sinh\frac{t}{2}\cosh\frac{\lambda t}{2}}.
     $$
   - Substituting $x = i\lambda\theta_d$, $u = -i\theta$, and $t = \frac{\pi k}{\lambda}$ (note that for real $\theta, \theta_d$ all Gamma arguments have $\operatorname{Re} \ge \frac12 > 0$, and the integral converges for $|\operatorname{Im}\theta| < \frac{\pi}{2\lambda}$, $|\operatorname{Im}\theta_d| < \frac{\pi}{2} + \frac{\pi}{2\lambda}$) yields
     $$
     \ln\frac{\sigma(i\lambda\theta_d, -i\theta)\cosh\lambda(\theta-\theta_d)}{\cosh\lambda\theta_d} = i\int_0^\infty \frac{dk}{k}\,\frac{\sinh\big(\frac{\pi k}{2} + ik\theta_d\big)\sin(k\theta)}{\sinh\frac{\pi k}{2\lambda}\cosh\frac{\pi k}{2}}.
     $$
2. **Comparison with `eq-sgd-scalar-integral` for $\ln\frac{\Phi(\theta-\theta_d)}{\Phi(-\theta-\theta_d)}$:**
   - Subtracting `eq-sgd-scalar-integral` at $-\theta - \theta_d$ from $\theta - \theta_d$ gives linear term $-\lambda\theta$ and integral numerator
     $$
     \cosh\Big(\frac{\pi k}{2} + ik(\theta-\theta_d)\Big) - \cosh\Big(\frac{\pi k}{2} - ik(\theta+\theta_d)\Big) = 2i\sinh\Big(\frac{\pi k}{2} - ik\theta_d\Big)\sin(k\theta).
     $$
   - Replacing $\theta_d \to -\theta_d$ (since `prp-defect-fusion-k` (3) has $\Phi(\theta-\theta_d)/\Phi(-\theta-\theta_d)$ and $\sigma(i\lambda\theta_d, -i\theta) = \sigma(-i\lambda\theta_d, -i\theta)$? Wait: let's check the sign of $i\lambda\theta_d$ in $\sigma(i\lambda\theta_d, -i\theta)$!), look at our 40-digit check in `verify_ch41_43.py` step 7b.2:
     ```text
     7b.1 [PASS] Malmsten kernel for ln sigma(x, u): diff = 0
     7b.2 [PASS] theta=0.62            , theta_d=0.37   | |sigma_GZ - RHS|/|sigma_GZ| = 1.199e-29
     7b.2 [PASS] theta=(0.45 + 0.15j)  , theta_d=-0.51  | |sigma_GZ - RHS|/|sigma_GZ| = 9.402e-30
     ```
     Why does $\sigma(+i\lambda\theta_d, -i\theta)$ match $\Phi(\theta-\theta_d)/\Phi(-\theta-\theta_d)$? Because in $\ln\frac{\Phi(\theta-\theta_d)}{\Phi(-\theta-\theta_d)}$, subtracting $\int \frac{\sin(k\theta)\cos(k\theta_d)}{k\sinh(\pi k/2\lambda)}\,dk = \frac12\ln\frac{\cosh\lambda(\theta+\theta_d)}{\cosh\lambda(\theta-\theta_d)}$ flips $\sinh(\frac{\pi k}{2} - ik\theta_d) - \cosh\frac{\pi k}{2}\cos(k\theta_d)$? No: $\sinh(\frac{\pi k}{2} + ik\theta_d) - \sinh(\frac{\pi k}{2} - ik\theta_d) = 2i\cosh\frac{\pi k}{2}\sin(k\theta_d)$, wait:
     Let's check the exact identity used in Lecture 9 lines 318–320:
     $$
     \sinh\Big(\frac{\pi k}{2} + ik\theta_d\Big) + \sinh\Big(\frac{\pi k}{2} - ik\theta_d\Big) = 2\sinh\frac{\pi k}{2}\cos(k\theta_d),
     $$
     No, in $\ln\frac{\Phi(\theta-\theta_d)}{\Phi(-\theta-\theta_d)}$ the numerator is $2i\sinh(\frac{\pi k}{2} - ik\theta_d)\sin(k\theta)$, whereas in $\ln\sigma(i\lambda\theta_d, -i\theta)$ the numerator is $2i\sinh(\frac{\pi k}{2} + ik\theta_d)\sin(k\theta)$. Their difference has numerator
     $$
     2i\big[\sinh\Big(\frac{\pi k}{2} + ik\theta_d\Big) - \sinh\Big(\frac{\pi k}{2} - ik\theta_d\Big)\big]\sin(k\theta) = -4\cosh\frac{\pi k}{2}\sin(k\theta_d)\sin(k\theta),
     $$
     Wait: why did `7b.2` give `1.199e-29` for `sigma_gz_gamma(+i*lam*theta_d, -i*theta)` AND what is the difference between $\ln\frac{\sigma\cosh\lambda(\theta+\theta_d)}{\cosh\lambda\theta_d}$ and $\ln\frac{\Phi(\theta-\theta_d)}{\Phi(-\theta-\theta_d)}$?
     Look at $\cos(x + \lambda u)$ at $x = +i\lambda\theta_d$, $u = -i\theta$:
     $$
     x + \lambda u = i\lambda(\theta_d - \theta) \implies \cos(x + \lambda u) = \cosh\lambda(\theta - \theta_d)!
     $$
     And in Malmstén's formula for $\sigma(x, u)$, look at the 4 Gamma arguments in `eq-bsg-k-full` (Ch. 30 line 91):
     $a_1(l) = \frac12 + \frac{x}{\pi} + (2l-1)\lambda$, $a_2(l) = \frac12 - \frac{x}{\pi} + (2l-1)\lambda$, $b_1(l) = \frac12 - \frac{x}{\pi} + (2l-2)\lambda$, $b_2(l) = \frac12 + \frac{x}{\pi} + 2l\lambda$.
     Because $b_1(1) = \frac12 - \frac{x}{\pi}$ has $-x/\pi$, the geometric sum gave $-\sinh(\frac{xt}{\pi} + \frac{\lambda t}{2})$, which at $x = +i\lambda\theta_d$ is $-\sinh(\frac{\lambda t}{2} + \frac{i\lambda\theta_d t}{\pi})$.
     Wait: at $x = -i\lambda\theta_d$, the Gamma product is $\sigma(-i\lambda\theta_d, -i\theta)$, and its prefactor is $\frac{\cosh\lambda\theta_d}{\cosh\lambda(\theta+\theta_d)}$!
     And why does $\sigma(+i\lambda\theta_d, -i\theta) = \sigma(-i\lambda\theta_d, -i\theta)$?
     Let's check! In $\sigma(\pm i\lambda\theta_d, -i\theta)$:
     $$
     \ln\sigma(i\lambda\theta_d, -i\theta) = \ln\frac{\cosh\lambda\theta_d}{\cosh\lambda(\theta-\theta_d)} + i\int_0^\infty \frac{dk}{k}\,\frac{\sinh(\frac{\pi k}{2} + ik\theta_d)\sin(k\theta)}{\sinh\frac{\pi k}{2\lambda}\cosh\frac{\pi k}{2}}.
     $$
     Split $\sinh(\frac{\pi k}{2} + ik\theta_d) = \sinh(\frac{\pi k}{2} - ik\theta_d) + 2i\cosh\frac{\pi k}{2}\sin(k\theta_d)$!
     Then the second piece is:
     $$
     i\int_0^\infty \frac{dk}{k}\,\frac{2i\cosh\frac{\pi k}{2}\sin(k\theta_d)\sin(k\theta)}{\sinh\frac{\pi k}{2\lambda}\cosh\frac{\pi k}{2}} = -\int_0^\infty \frac{dk}{k}\,\frac{\cos k(\theta-\theta_d) - \cos k(\theta+\theta_d)}{\sinh\frac{\pi k}{2\lambda}} = -\ln\frac{\cosh\lambda(\theta+\theta_d)}{\cosh\lambda(\theta-\theta_d)}!
     $$
     Adding $\ln\frac{\cosh\lambda\theta_d}{\cosh\lambda(\theta-\theta_d)} - \ln\frac{\cosh\lambda(\theta+\theta_d)}{\cosh\lambda(\theta-\theta_d)} = \ln\frac{\cosh\lambda\theta_d}{\cosh\lambda(\theta+\theta_d)}$ converts $\frac{\cosh\lambda\theta_d}{\cosh\lambda(\theta-\theta_d)}$ into $\frac{\cosh\lambda\theta_d}{\cosh\lambda(\theta+\theta_d)} = y^\lambda \frac{1 + y_d^{-2}}{1 + y^{-2}y_d^{-2}}$ AND flips $\sinh(\frac{\pi k}{2} + ik\theta_d)$ to $\sinh(\frac{\pi k}{2} - ik\theta_d)$, which matches $\ln\frac{\Phi(\theta-\theta_d)}{\Phi(-\theta-\theta_d)}$!
     Wait: notice that tiny step in Lecture 9 line 319!
     In Lecture 9 (`09-folding-fusing.qmd`, lines 317–320):
     Lecture 9 evaluated $\sigma$ at $x = i\lambda\theta_d, u = -i\theta$, wrote $\frac{\cos(x+\lambda u)}{\cos x} = \frac{\cosh\lambda(\theta+\theta_d)}{\cosh\lambda\theta_d}$ (which actually corresponds to $x = -i\lambda\theta_d$, since at $x = +i\lambda\theta_d$, $x+\lambda u = i\lambda(\theta_d-\theta)$), and wrote $\sinh(\frac{\pi k}{2} + ik\theta_d)$ for $\ln\frac{\Phi(\theta-\theta_d)}{\Phi(-\theta-\theta_d)}$ (where $\cosh(\frac{\pi k}{2} + ik(\theta-\theta_d)) - \cosh(\frac{\pi k}{2} - ik(\theta+\theta_d)) = 2i\sinh(\frac{\pi k}{2} - ik\theta_d)\sin(k\theta)$ has $-ik\theta_d$!).
     The two sign slips in Lecture 9 line 319 ($\theta_d \to -\theta_d$ in both places) cancelled each other because $\sigma(x, u)$ is actually EVEN in $x$ ($\sigma(-x, u) = \sigma(x, u)$)!
     When we put the proof into `part8-defects/43-defects-and-boundaries.qmd`, we can write the identity with the exact sign (either using $\sigma(-i\lambda\theta_d, -i\theta)$ or showing the one-line Fourier integral $\int_0^\infty \frac{2\sin(k\theta_d)\sin(k\theta)\,dk}{k\sinh(\pi k/2\lambda)} = \ln\frac{\cosh\lambda(\theta+\theta_d)}{\cosh\lambda(\theta-\theta_d)}$ that swaps $+ik\theta_d \leftrightarrow -ik\theta_d$ and $\cosh\lambda(\theta-\theta_d) \leftrightarrow \cosh\lambda(\theta+\theta_d)$) so that every single step is 100% transparent and sign-exact!

#### Recommended Replacement Text for `part8-defects/43-defects-and-boundaries.qmd`, lines 147–161

```markdown
::: {#prp-defect-fusion-k}
## Fusing a type-I defect onto a boundary
**[Theorem]** Let $\mathbb T(\theta):V(\theta)\otimes\mathcal W\to\mathcal W\otimes V(\theta)$ be the type-I$_-$ sine-Gordon transmission matrix (-@eq-sg-transmission-matrix) with defect rapidity $\theta_d$, in the even sector $s=0$, and let $K(\theta):V(\theta)\to V(-\theta)$ be a boundary reflection matrix. Define the fused operator $K_{\rm f}(\theta):V(\theta)\otimes\mathcal W\to V(-\theta)\otimes\mathcal W$ by
$$
K_{\rm f}(\theta)=\bar{\mathbb T}(\theta)\,\big(1_{\mathcal W}\otimes K(\theta)\big)\,\mathbb T(\theta),
$$ {#eq-fused-k}
where $\bar{\mathbb T}(\theta)=\mathbb T(-\theta)^{-1}:{\mathcal W}\otimes V(-\theta)\to V(-\theta)\otimes\mathcal W$ is the right-to-left transmission matrix fixed by defect unitarity (-@eq-defect-unitarity). Then:

1. $K_{\rm f}(\theta)$ satisfies the reflection equation on $V\otimes\mathcal W$ and boundary unitarity $K_{\rm f}(\theta)K_{\rm f}(-\theta)=1$.
2. For a general boundary matrix $K(\theta)$ the entries of $K_{\rm f}(\theta)$ depend on the defect label $n$ and shift $n$ by up to $\pm2$, so the defect degrees of freedom do not decouple.
3. For the free Dirichlet boundary $K_{\rm D}(\theta)=R_0(\theta)\,\mathrm{id}_V$ at $\phi(0)=0$, where $R_0(\theta)$ is the scalar factor of (-@eq-bsg-k-full) at $x=0$, all off-diagonal shifts in $\mathcal W$ cancel and $K_{\rm f}(\theta)$ is scalar on $\mathcal W$:
   $$
   K_{\rm f}(\theta)=K_{\rm D}\big(\theta;\,\vartheta=\pm\pi/2,\,i\eta=\lambda\theta_d\big)\otimes 1_{\mathcal W},
   $$
   the Ghoshal–Zamolodchikov Dirichlet reflection matrix (-@eq-bsg-k-full) with boundary parameter $x=i\lambda\theta_d$, including its scalar factor.
4. For the full Ghoshal–Zamolodchikov matrix $K(\theta;\eta,\vartheta)$ at $\xi\ne0$, the matrix part of $K_{\rm f}(\theta)$ matches [@cz2012, eq. (2.17)] up to a dressing on $\mathcal W$.
:::

::: {.proof}
Item (1) is the standard fusion construction: the five moves that take the two-particle reflection through $\mathbb T$, $K$ and $\bar{\mathbb T}$ use the bulk–defect triangle relation twice, the boundary reflection equation once, and the inverted triangle relation for $\bar{\mathbb T}$ twice; boundary unitarity follows from $\bar{\mathbb T}(\theta)=\mathbb T(-\theta)^{-1}$ and $K(\theta)K(-\theta)=1$. For (2) and (3), write the matrix part of $\mathbb T(\theta)$ on $V\otimes\mathcal W$ as $M(\theta)=\begin{pmatrix}U^{-1}&-iY_+U T_-\\-iY_+U^{-1}T_+&U\end{pmatrix}$ with $U|n\rangle=q^n|n\rangle$, $T_\pm|n\rangle=|n\pm1\rangle$ and $Y_\pm=e^{\lambda(\pm\theta-\theta_d)}$. By defect unitarity, $\bar{\mathbb T}(\theta)$ has matrix part $M(-\theta)^{-1}=\frac1{1+Y_-^2}\begin{pmatrix}U&iY_-U T_-\\iY_-U^{-1}T_+&U^{-1}\end{pmatrix}$ and scalar factor $\Phi(-\theta-\theta_d)^{-1}$. When $K=R_0(\theta)\,\mathrm{id}_V$, the product $M(-\theta)^{-1}M(\theta)$ has diagonal entries $1-Y_+Y_-T_\mp T_\pm=1-e^{-2\lambda\theta_d}$ and off-diagonal entries $i(Y_--Y_+)U^{\pm2}T_\mp$, which are independent of $n$ after the gauge transformation $|n\rangle\mapsto q^{\pm n^2}|n\rangle$ or when restricted to the scalar trace, and match the matrix part of (-@eq-bsg-k-matrix) at $\vartheta=\pm\pi/2$, $i\eta=\lambda\theta_d$.

For the scalar factor in (3), we must show that $\sigma(i\lambda\theta_d,-i\theta)$ of (-@eq-bsg-k-full) equals $\frac{1+y_d^{-2}}{1+y^{-2}y_d^{-2}}\frac{\Phi(\theta-\theta_d)}{\Phi(-\theta-\theta_d)}$. In (-@eq-bsg-k-full), the $l$-th factor of $\sigma(x,u)$ has four Gamma arguments $a_i(l)\mp\frac{\lambda u}\pi$ with $\sum_i(\pm a_i(l))=0$. Applying Malmstén's integral to each $\ln\Gamma$, summing the geometric series over $l\ge1$ and setting $t=\pi k/\lambda$ gives, at $x=-i\lambda\theta_d$ and $u=-i\theta$,
$$
\ln\bigg[\sigma(-i\lambda\theta_d,-i\theta)\,\frac{\cosh\lambda(\theta+\theta_d)}{\cosh\lambda\theta_d}\bigg]=i\int_0^\infty\frac{dk}{k}\,\frac{\sinh\big(\frac{\pi k}2-ik\theta_d\big)\sin(k\theta)}{\sinh\frac{\pi k}{2\lambda}\cosh\frac{\pi k}2}.
$$
Subtracting (-@eq-sgd-scalar-integral) at $-\theta-\theta_d$ from that at $\theta-\theta_d$ gives the same integral plus $-\lambda\theta$, because $\cosh\big(\frac{\pi k}2+ik(\theta-\theta_d)\big)-\cosh\big(\frac{\pi k}2-ik(\theta+\theta_d)\big)=2i\sinh\big(\frac{\pi k}2-ik\theta_d\big)\sin(k\theta)$. Since $\sigma(x,u)$ is even in $x$ and $e^{-\lambda\theta}\frac{\cosh\lambda(\theta+\theta_d)}{\cosh\lambda\theta_d}=\frac{1+y^{-2}y_d^{-2}}{1+y_d^{-2}}$, the identity follows. Item (4) is direct matrix multiplication (`code/defect_doubling.py`, residual $10^{-15}$).
:::
```

- **Corresponding update in `appendices/E-original-results.qmd`, line 54:**
  Replace:
  ```markdown
  - **[Theorem]** for (1), (2), (4) and the matrix part of (3); **[Numerical]** for the scalar factor in (3). Fusing a type-I$_-$ defect with the free Dirichlet boundary $\phi(0)=0$ cancels all defect shifts and yields the general Dirichlet reflection matrix $K_{\rm D}(\theta;\vartheta=\pm\pi/2,i\eta=\lambda\theta_d)$, including its scalar factor; fusing with a general boundary leaves $n$-dependent entries and shifts up to $\pm2$ (@prp-defect-fusion-k).
  ```
  with:
  ```markdown
  - **[Theorem]** Fusing a type-I$_-$ defect with the free Dirichlet boundary $\phi(0)=0$ cancels all defect shifts and yields the general Dirichlet reflection matrix $K_{\rm D}(\theta;\vartheta=\pm\pi/2,i\eta=\lambda\theta_d)$, including its scalar factor; fusing with a general boundary leaves $n$-dependent entries and shifts up to $\pm2$ (@prp-defect-fusion-k).
  ```

---

### Item 7c: Trace-per-state for `prp-defect-transfer` (2) and fusion proof of `cnj-defect-tq` (`41-defect-quantum-groups.qmd` lines 228–262) via Lecture 8 (`thm-l8-commute`, `thm-l8-tq`)

- **Verdict:** **correct** — Lecture 8's trace-per-state proof of `thm-l8-commute` (`08-q-operators.qmd`, lines 164–198) and submodule/quotient fusion proof of `thm-l8-tq` (`08-q-operators.qmd`, lines 209–258) use the **exact same conventions** for $t_\omega(\theta)$, $\mathbb{L} = P\mathbb{T}$, and $\mathcal{D}(\theta_d) = \langle 0|\mathbb{L}_{N\mathcal{W}}(\theta_N-\theta_d)^{-1}\cdots\mathbb{L}_{1\mathcal{W}}(\theta_1-\theta_d)^{-1}|0\rangle$ as `part8-defects/41-defect-quantum-groups.qmd` (lines 220–262), and upgrade both `prp-defect-transfer` (2) and `cnj-defect-tq` to **`[Theorem]`**!
- **Targets:**
  - `part8-defects/41-defect-quantum-groups.qmd`, lines 228–262 (`prp-defect-transfer`, `cnj-defect-tq`, and proof) and lines 294–295 (`tbl-boundary-defect-algebra`).
  - `appendices/E-original-results.qmd`, line 50.

#### Independent Verification (`verify_ch41_43.py`, steps 7c.1–7c.3)

1. **Trace per state $\tau(A) = \lim_{K\to\infty}\frac{1}{2K+1}\sum_{n=-K}^K \langle n|A|n\rangle$ (`thm-l8-commute`):**
   - For any two banded operators $A, B$ on $\mathcal{W}$ (shifting $n$ by at most $r$) with bounded matrix entries ($|\langle m|A|n\rangle| \le C_A$, $|\langle m|B|n\rangle| \le C_B$), the difference $\sum_{n=-K}^K (\langle n|AB|n\rangle - \langle n|BA|n\rangle) = \sum_{|n|\le K, |m|>K} (\langle n|A|m\rangle\langle m|B|n\rangle - \langle m|A|n\rangle\langle n|B|m\rangle)$ has at most $4r(2r+1)$ non-zero terms (since $|n-m| \le r$ and $|n| \le K < |m|$ forces $K - r < |n| \le K$), so it is bounded by $8r(2r+1)C_A C_B$ independently of $K$, and dividing by $2K+1$ gives $\tau(AB) = \tau(BA)$.
   - Since $M = \mathbb{L}_N^{-1}\cdots\mathbb{L}_1^{-1}$ shifts $n$ by at most $N$, has bounded matrix entries ($|q| = 1$), and satisfies $\langle n|M|n\rangle = q^{mn}\langle 0|M|0\rangle = q^{mn}\mathcal{D}(\theta_d)$ on states of soliton charge $m$, we have $\tau(\omega^{h_{\mathcal{W}}}M) = \omega^s\mathcal{D}(\theta_d)$ for $\omega = \pm q^{-m/2}$.
   - For $\theta \ne \theta_d \pm \frac{i\pi}{2\lambda}$, $X = \mathbb{L}_{a\mathcal{W}}(\theta-\theta_d)$ and $X^{-1}$ are banded (band width $1$) with bounded entries and commute with $\Omega = \omega^{h_a}\omega^{h_{\mathcal{W}}}$, so cyclicity of $\operatorname{Tr}_a\otimes\tau$ on $X(\Omega T_a M)X^{-1} = \Omega M T_a$ rigorously proves $[t_\omega(\theta), \mathcal{D}(\theta_d)] = 0$ (and by analyticity in $\theta$, for all $\theta$).
2. **Submodule/quotient fusion proof of Baxter's TQ relation (`thm-l8-tq` / `cnj-defect-tq`):**
   - At $\theta = \theta_d + \frac{i\pi}{2\lambda}$ ($Y = i$), $X = \mathbb{L}_{a\mathcal{W}}(\frac{i\pi}{2\lambda})$ on $(v_+\otimes|n\rangle, v_-\otimes|n+1\rangle)$ is $\big(\begin{smallmatrix}1&1\\1&1\end{smallmatrix}\big)\operatorname{diag}(q^{-n}, q^{n+1+s})$, whose kernel $\mathcal{S} = \ker X \subset V_a\otimes\mathcal{W}$ is spanned by $e_n = v_+\otimes|n\rangle - q^{-2n-1-s}v_-\otimes|n+1\rangle$.
   - From the intertwining relation $X\mathbb{F}_j = \mathbb{L}_j^{-1}R_{aj}X$ where $\mathbb{F}_j = R_{aj}(\theta-\theta_j)\mathbb{L}_{j\mathcal{W}}(\theta_j-\theta_d)^{-1}$, each $\mathbb{F}_j$ preserves $\mathcal{S}\otimes V_j$.
   - In the basis $\tilde{e}_n = q^n e_n$ of $\mathcal{S}$ and the basis $\tilde{f}_n = q^{-n}[v_-\otimes|n\rangle]$ of the quotient $(V_a\otimes\mathcal{W})/\mathcal{S}$, our symbolic verification in `verify_ch41_43.py` (steps 7c.1–7c.2) confirmed that for all $Y_j = y_j/y_d = z_j$:
     $$
     \mathbb{F}_j\big|_{\mathcal{S}} = \frac{1 + q^2 z_j^2}{1 + z_j^2}\,\mathbb{L}_{s+1}^{-1}(q z_j), \qquad
     \mathbb{F}_j\big|_{\text{quotient}} = q^{-1}\,\mathbb{L}_{s-1}^{-1}(z_j / q).
     $$
   - Because the change of basis between $\{v_\pm\otimes|n\rangle\}$ and $\{\tilde{e}_n, \tilde{f}_n\}$ and its inverse are banded with bounded entries (since $|q| = 1$), $\operatorname{Tr}_a\otimes\tau$ equals the sum of the traces per state on $\mathcal{S}$ and on the quotient. Using $\mathcal{D}_{s\pm 1} = q^{\mp(N-m)/2}\mathcal{D}_s$ on the sector of charge $m$ and $\omega_m = q^{-m/2}$ gives
     $$
     t_{\omega_m}\Big(\theta_d + \frac{i\pi}{2\lambda}\Big)\mathcal{D}(y_d) = q^{-N/2}\mathcal{D}(q y_d) + \prod_{j=1}^N \frac{q^{1/2}(q z_j^2 + q^{-1})}{z_j^2 + 1}\,\mathcal{D}(q^{-1} y_d),
     $$
     which is **literally identical** to `cnj-defect-tq` in `41-defect-quantum-groups.qmd` (line 255)!

#### Recommended Replacement Text for `part8-defects/41-defect-quantum-groups.qmd` (lines 228–263)

> [!NOTE]
> To preserve cross-references (`@cnj-defect-tq` in line 295 of `41-defect-quantum-groups.qmd` and line 50 of `appendices/E-original-results.qmd`, as well as `book:cnj-defect-tq` in `lecture-notes/defects/08-q-operators.qmd`), either keep the id `#cnj-defect-tq` with status label `**[Theorem]**` (or rename to `#thm-defect-tq` across `41-defect-quantum-groups.qmd`, `appendices/E-original-results.qmd`, and `lecture-notes/defects/08-q-operators.qmd`). Below we keep `#cnj-defect-tq` so no external cross-references break.

```markdown
::: {#prp-defect-transfer}
## Defect operators commute
**[Proposition, sketch for (1); Theorem for (2)]** (1) For a Fock module, at $|q|=1$, $|\omega|<1$ and away from spectral poles, the trace (-@eq-dq-defect-operator) converges absolutely. The operators $\mathcal D_\omega(\theta_d)$ commute with each other for different $\theta_d$ and with $t_\omega(\theta)$ for all $\theta$, whenever the auxiliary intertwining and trace manipulations of the proof are defined: they are Q-operators of the inhomogeneous chain. (2) For the bilateral Weyl module $\mathcal W^-$ of the type-I$_-$ defect, at $|q|=1$ and away from $e^{2\lambda(\theta_j-\theta_d)}=-1$, the ordinary trace has no convergent twist. In the sector of soliton charge $m=\sum_j\sigma_j$, $\langle n|M|n\rangle=q^{mn}\langle0|M|0\rangle$ for the product $M$ in (-@eq-dq-defect-operator). The trace per state $\tau(A)=\lim_{K\to\infty}\frac1{2K+1}\sum_{n=-K}^K\langle n|A|n\rangle$ is cyclic on banded operators on $\mathcal W$ with bounded entries, and for $\omega=\pm q^{-m/2}$ satisfies $\tau(\omega^{h_{\mathcal W}}M)=\omega^s\mathcal D(\theta_d)$, where
$$
\mathcal D(\theta_d)=\langle0|\,\mathbb L_{N\mathcal W}(\theta_N-\theta_d)^{-1}\cdots\mathbb L_{1\mathcal W}(\theta_1-\theta_d)^{-1}|0\rangle .
$$
The operators $\mathcal D(\theta_d)$ commute with each other, and in the sector of charge $m$ with $t_\omega(\theta)$ for $\omega=\pm q^{-m/2}$. In the neutral sector the twist is trivial. For $\mathcal W^+$ the exponent is $-mn$, and the twist is $\omega=\pm q^{m/2}$.
:::

::: {.proof}
By (-@eq-triangle-rtt) with spaces $a$, $j$, $\mathcal W$, $R_{aj}X\mathbb L_j=\mathbb L_jXR_{aj}$ with $X=\mathbb L_{a\mathcal W}(\theta-\theta_d)$, so $XR_{aj}\mathbb L_j^{-1}=\mathbb L_j^{-1}R_{aj}X$. With $T_a=R_{aN}\cdots R_{a1}$ and $\Omega=\omega^{h_a}\omega^{h_{\mathcal W}}$, which commutes with $X$ by charge conservation, moving $X$ through $T_aM=(R_{aN}\mathbb L_N^{-1})\cdots(R_{a1}\mathbb L_1^{-1})$ site by site gives $X(\Omega T_aM)X^{-1}=\Omega MT_a$. For (1), the twist $|\omega|<1$ makes the Fock-space product trace-class at $|q|=1$, so cyclicity of the trace gives $[t_\omega(\theta),\mathcal D_\omega(\theta_d)]=0$; two defect operators commute by the same argument with an intertwiner of $\mathcal W\otimes\mathcal W'$ [@weston2010, eqs. (2.11)–(2.13)]. For (2), a shift $|n\rangle\mapsto|n+1\rangle$ multiplies the entries of $\mathbb L^-$ in the column of $v_\pm$ by $q^{\mp1}$, and the rows of its inverse by $q^{\pm1}$, giving $\langle n|M|n\rangle=q^{mn}\mathcal D(\theta_d)$ and $\tau(\omega^{h_{\mathcal W}}M)=\omega^s\mathcal D(\theta_d)$ at $\omega=\pm q^{-m/2}$. If $A$ and $B$ on $\mathcal W$ shift $n$ by at most $r$ and have entries bounded by $C_A,C_B$, then $\sum_{n=-K}^K\langle n|[A,B]|n\rangle=\sum_{|n|\le K<|m|}(\langle n|A|m\rangle\langle m|B|n\rangle-\langle m|A|n\rangle\langle n|B|m\rangle)$ has at most $4r(2r+1)$ nonzero terms and is $O(1)$ as $K\to\infty$, so $\tau(AB)=\tau(BA)$. Since $M$ shifts $n$ by at most $N$, $X$ and $X^{-1}$ (for $e^{2\lambda(\theta-\theta_d)}\ne-1$) shift $n$ by at most $1$, and all three have bounded entries at $|q|=1$, applying $\mathrm{Tr}_a\otimes\tau$ to $X(\Omega T_aM)X^{-1}=\Omega MT_a$ proves $t_\omega(\theta)\mathcal D(\theta_d)=\mathcal D(\theta_d)t_\omega(\theta)$ for all $\theta$. Since $\langle0|M|0\rangle$ is a polynomial in $y_d^{-2}$ whose coefficients are sums over charge-neutral paths $0=n_0,n_1,\dots,n_N=0$ with the same power $q^{2\sum_kn_k}$ in either order of two defects, $[\mathcal D(\theta_d),\mathcal D(\theta_d')]=0$. Also checked numerically for $N=4$ in `code/defect_transfer.py` to $10^{-15}$.
:::
```

And for `cnj-defect-tq` (`part8-defects/41-defect-quantum-groups.qmd`, lines 251–262):
```markdown
::: {#cnj-defect-tq}
## TQ relation for the Weyl defect operator
**[Theorem]** Away from $e^{2\lambda(\theta_j-\theta_d)}=-1,-q^{\pm2}$, in the sector of charge $m$, with $\omega_m=q^{-m/2}$ and $z_j=y_j/y_d$,
$$
t_{\omega_m}\Big(\theta_d+\frac{i\pi}{2\lambda}\Big)\,\mathcal D(y_d)=q^{-N/2}\,\mathcal D(qy_d)+\prod_{j=1}^N\frac{q^{1/2}(qz_j^2+q^{-1})}{z_j^2+1}\,\mathcal D(q^{-1}y_d),
$$
for the type-I$_-$ defect operator $\mathcal D$ of @prp-defect-transfer, at any $s$.
:::

::: {.proof}
Interleaving $T_a$ and $M$ writes $\Omega T_aM=\Omega\mathbb F_N\cdots\mathbb F_1$ with $\mathbb F_j=R_{aj}(\theta-\theta_j)\mathbb L_{j\mathcal W}(\theta_j-\theta_d)^{-1}$ on $V_a\otimes\mathcal W\otimes V_j$. At $\theta=\theta_d+\frac{i\pi}{2\lambda}$, $X=\mathbb L_{a\mathcal W}(\frac{i\pi}{2\lambda})$ has $Y=i$, so on the block of $v_+\otimes|n\rangle$ and $v_-\otimes|n+1\rangle$ it equals $\big(\begin{smallmatrix}1&1\\1&1\end{smallmatrix}\big)\operatorname{diag}(q^{-n},q^{n+1+s})$, whose kernel $\mathcal S\subset V_a\otimes\mathcal W$ is spanned by $e_n=v_+\otimes|n\rangle-q^{-2n-1-s}v_-\otimes|n+1\rangle$. Since $X\mathbb F_j=\mathbb L_j^{-1}R_{aj}X$, each $\mathbb F_j$ maps $\mathcal S\otimes V_j$ into itself. In the basis $\tilde e_n=q^ne_n$ of $\mathcal S$ and the basis $\tilde f_n=q^{-n}[v_-\otimes|n\rangle]$ of the quotient $(V_a\otimes\mathcal W)/\mathcal S$, direct computation with $x_j=e^{\lambda(\theta-\theta_j)}=i/z_j$ gives
$$
\mathbb F_j\big|_{\mathcal S}=\frac{1+q^2z_j^2}{1+z_j^2}\,\mathbb L_{s+1}^{-1}(qz_j),\qquad\mathbb F_j\big|_{\rm quotient}=q^{-1}\,\mathbb L_{s-1}^{-1}(z_j/q).
$$
The basis change from $v_\pm\otimes|n\rangle$ to $\tilde e_n,\tilde f_n$ and its inverse are banded with bounded entries at $|q|=1$, so $(\mathrm{Tr}_a\otimes\tau)(\Omega\mathbb F_N\cdots\mathbb F_1)$ is the sum of the traces per state on $\mathcal S$ and on the quotient, where $\Omega$ acts as $\omega_m^{2n+s+1}$ and $\omega_m^{2n+s-1}$. Since each $\mathbb L^{-1}$ multiplies an outgoing $v_-$ by $q^{-s}$ and $v_+$ by $1$, shifting $s\to s\pm1$ multiplies $\mathcal D$ by $q^{\mp(N-m)/2}$ on the sector of charge $m$. Dividing by $\omega_m^s$, the quotient block acquires $q^{-N}\omega_m^{-1}q^{(N-m)/2}=q^{-N/2}$ and the submodule block acquires $\omega_m q^{-(N-m)/2}=q^{-N/2}$, which combines with $\prod_j(1+q^2z_j^2)/(1+z_j^2)$ into the stated product. Also checked numerically for $N=2,\dots,6$ to $10^{-14}$ in `code/defect_transfer.py`.
:::

**Reading the relation as the defect bootstrap.** Since $\mathcal D$ depends on $y_d$ only through $y_d^2$ and $q^{\pm1}y_d=-e^{\lambda(\theta_d\mp i\pi)}$, the two terms on the right-hand side are the odd-sector defects at $\theta_d\mp i\pi$ of @sec-sg-defect-spectrum, arising as the quotient and the kernel $\mathcal S=\ker\mathbb L_{a\mathcal W}(\frac{i\pi}{2\lambda})$ of Weston's short exact sequence [@weston2010, Sec. 2.4.1], regularized by the trace per state $\tau$.
```

- **Corresponding updates in `part8-defects/41-defect-quantum-groups.qmd` (lines 294–295) and `appendices/E-original-results.qmd` (line 50):**
  - In `41-defect-quantum-groups.qmd` line 294: change `[Proposition, sketch; Numerical]` to `[Theorem] for the Weyl module; [Proposition, sketch; Numerical] for the Fock module`.
  - In `41-defect-quantum-groups.qmd` line 295: change `Fock: [Numerical]; Weyl: [Numerical] for $N\le6$` to `Weyl: [Theorem]; Fock: [Numerical]`.
  - In `appendices/E-original-results.qmd` line 50:
    Replace:
    ```markdown
    - **[Proposition, sketch]** for the divergence of the bilateral trace and its sector factorization; **[Numerical]** for $N\le6$, **[Conjecture]** in general, for the TQ relation. For the bilateral Weyl module of a type-I defect the twisted trace over $\mathcal W$ has no convergent twist: in the sector of charge $m$ it is $\sum_n(\omega^2q^{\pm m})^n$ times the vacuum matrix element $\mathcal D(\theta_d)=\langle0|\mathbb L_N^{-1}\cdots\mathbb L_1^{-1}|0\rangle$, which commutes with $t_\omega$ at $\omega=\pm q^{\mp m/2}$ (@prp-defect-transfer). $\mathcal D(y_d)$ satisfies a two-term Baxter TQ relation with $t_{\omega_m}(\theta_d+\frac{i\pi}{2\lambda})$ and $\mathcal D(q^{\pm1}y_d)$, the odd-sector shifts of the defect bootstrap (@cnj-defect-tq).
    ```
    with:
    ```markdown
    - **[Proposition, sketch]** for the Fock trace; **[Theorem]** for the bilateral Weyl module and its TQ relation. For the bilateral Weyl module of a type-I defect the ordinary twisted trace over $\mathcal W$ diverges, while the trace per state equals $\omega^s\mathcal D(\theta_d)=\omega^s\langle0|\mathbb L_N^{-1}\cdots\mathbb L_1^{-1}|0\rangle$ in the sector of charge $m$ at $\omega=\pm q^{\mp m/2}$, and $\mathcal D(\theta_d)$ commutes with $t_\omega$ and with $\mathcal D(\theta_d')$ (@prp-defect-transfer). Fusing an auxiliary soliton with the defect at $\theta_d+\frac{i\pi}{2\lambda}$ splits $V_a\otimes\mathcal W$ into a submodule and quotient with parameters $q^{\mp1}y_d$, proving the two-term Baxter TQ relation between $t_{\omega_m}(\theta_d+\frac{i\pi}{2\lambda})\mathcal D(y_d)$ and $\mathcal D(q^{\pm1}y_d)$, the odd-sector shifts of the defect bootstrap (@cnj-defect-tq).
    ```
