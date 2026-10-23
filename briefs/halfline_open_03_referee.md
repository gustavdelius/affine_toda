# Referee report: half-line open problem 03
Date: 2026-10-11. Referee: Gemini 3.8 Flash (High), via Antigravity. Tier: B.

## Summary
The worker's report succeeds in deriving the Neumann particle reflection amplitudes of $a_n^{(1)}$ from the $s=+1$ non-diagonal soliton reflection matrix. Every step of the derivation—the boundary unitarity condition $K(-\theta)^{\rm T}K(\theta) = U \mathbf 1$, the crossing function $\Phi_+$, the $B_1$ singlet projection $c_+$, the minimal scalar factor $A^{(+)}$, the CDD factor, and the continuation $\omega \to -2/B$ to $K^N_1$—has been independently verified. The report corrects a significant premise of the brief and the textbook regarding boundary unitarity for non-diagonal reflection matrices, resolves the branch convention of the coideal matrix, and provides an exact analytic proof of the unitarity function $U$. All numerical checks reproduce the claimed closed forms to machine precision or the stated product truncation precision, and the continued amplitudes satisfy all scattering axioms and match Kim's one-loop results. The open problem is decisively advanced to [Numerical] status for $n \le 5$.

## Verdicts
| # | Claim (short) | Label claimed | Verdict | Label you recommend |
|---|---|---|---|---|
| 1 | Boundary unitarity is $A(\theta)A(-\theta)K(-\theta)^{\rm T}K(\theta)=1$, not $K(\theta)K(-\theta)$ | [Correction] | CONFIRMED | [Theorem] |
| 2 | Branch in book's `k_closed` selects wrong sign class for odd $N$ and $\cos(\pi\omega/2)<0$ | [Correction] | CONFIRMED | [Correction] |
| 3 | Unitarity function $K(-\theta)^{\rm T}K(\theta)=U\mathbf 1$ with $U=(s\cosh2T\theta-\cos\pi T)/(2\sin^2\frac{\pi\omega}2)$ | [Theorem] | CONFIRMED | [Theorem] |
| 4 | Crossing function $\Phi_+=\sin\pi(\frac\omega2-\mu)/\sin\pi(\frac\omega2+\mu)$, $\Phi_-=\cos\pi(\mu-\frac\omega2)/\cos\pi(\mu+\frac\omega2)$ | [Numerical] | CONFIRMED | [Numerical] ($n\le5$), [Theorem] ($n=2$) |
| 5 | Singlet projection $(1\otimes K)\check R(1\otimes K^{\rm T})$ maps $B_1$ singlet to $c_\pm(\mu)$ times itself | [Numerical] | CONFIRMED | [Numerical] |
| 6 | Scalar factors $A^{(\pm)}$ solve unitarity and crossing; match GZ ($n=1$) and Gandenberger ($n=2$) | [Numerical] | CONFIRMED | [Numerical] |
| 7 | Soliton bootstrap $V^{\otimes n}\to\bar V$ holds for $A^{(+)}$ ($n\le5$) and $A^{(-)}$ ($n\le4$) | [Numerical] | CONFIRMED | [Numerical] |
| 8 | CDD factor analysis: $(\frac T2)_I^2$ fixes zero-momentum pole; $(\frac T2\mp\frac\omega2)_I$ forced for even $n$ | [Analysis] | CONFIRMED | [Analysis] |
| 9 | Main result: $K_{B_1}^{(+)}$ continues under $\omega\to-2/B$ to $K^N_1$, generating $K^N_a$ | [Numerical] | CONFIRMED | [Numerical] (solitons), [Theorem] (block identity) |
| 10 | The $s=-1$ bonus: $K_{B_1}^{(-)}$ continues to $K^{(-)}_a=-(a)(h-a)K_a^{\rm unif}$, matching Gandenberger (6.4) | [Numerical] | CONFIRMED | [Numerical] |
| 11 | Semiclassical time advance of $A^{(+)}$ matches Neumann classical limit with $1/\omega$ errors | [Numerical] | CONFIRMED | [Numerical] |
| 12 | Proposed replacement text for `sec-neumann-duality`, `sec-an-nondiagonal-k`, and open problem 3 | [Text] | CONFIRMED | [Text] |

## Details

### Claim 1
**Boundary unitarity for non-diagonal $K$.**
The reflection matrix $K: V(\theta) \to \bar V(-\theta)$ maps a particle in representation $V$ to an antiparticle in representation $\bar V$. A second reflection off the boundary is therefore governed by the reflection matrix of the antiparticle, $K_{\bar V}(-\theta)$.
Under charge conjugation $J$ (the diagram automorphism $j \leftrightarrow n+1-j$ which fixes node 0 and preserves the boundary coideal parameters $\hat\epsilon_j = \epsilon_*$), $K_{\bar V} = J K J = K^{\rm T}$.
Applying the two reflections successively gives $K_{\bar V}(-\theta) K(\theta) = K(-\theta)^{\rm T} K(\theta)$.
We verified numerically in Python that $K(-\theta) K(\theta)$ is not scalar for $n \ge 2$ (e.g. at $n=2$, off-diagonal elements are nonzero, $\approx -0.098 - 0.945j$), whereas $K(-\theta)^{\rm T} K(\theta)$ is strictly proportional to the identity matrix with off-diagonal elements zero to machine precision ($< 10^{-15}$).
The physical unitarity condition is therefore $A(\theta)A(-\theta) K(-\theta)^{\rm T} K(\theta) = \mathbf 1$. Confirmed.

### Claim 2
**Branch in book code.**
In `code/boundary_coideal.py`, `k_closed` computes:
`w = np.sqrt(-qq + 0j)`
`es = 1/np.sqrt((1 - qq)*(1 - 1/qq))`
At physical $q = -e^{-i\pi\omega}$, $-q = e^{-i\pi\omega}$. When $\omega > 1$ (e.g. $\omega = 1.83$), numpy's principal square root maps $\text{Arg}(-q) \in (-\pi, \pi]$ by halving the wrapped argument, which introduces a sign flip relative to the analytic branch $w = e^{-i\pi\omega/2}$. Furthermore, `es` takes a square root of a positive real number, yielding $\epsilon_* > 0$ regardless of $\omega$, whereas the analytic continuation requires $\epsilon_* = \frac{1}{2\cos(\pi\omega/2)}$, which changes sign when $\cos(\pi\omega/2) < 0$. The book scripts only test at $\omega = 0.37$, where this issue does not manifest. The worker's `k_book` in `nk_common.py` adopts the analytic definitions $w = e^{-i\pi\omega/2}$ and $\epsilon_* = \frac{1}{2\cos(\pi\omega/2)}$, which correctly preserves the sign class for all $\omega$. Confirmed.

### Claim 3
**The unitarity function $U_n$.**
The worker claims $K(-\theta)^{\rm T} K(\theta) = U \mathbf 1$ with $U = \frac{s\cosh 2T\theta - \cos\pi T}{2\sin^2(\pi\omega/2)}$.
We independently verified the complete analytic proof symbolically using SymPy:
Expressing $K$ in terms of the twisted shift operator $Z$ ($Z^N = s \mathbf 1$), the eigenvalue on an eigenvector of $Z$ with eigenvalue $\zeta$ ($\zeta^N = s$) is:
$$\kappa(\zeta) = \frac{(1 - s X^N)(X\zeta + q)}{w X^{N/2}(w - w^{-1})(X\zeta - 1)} .$$
For $K(-\theta)^{\rm T}$, $X \to w^2/X$ and $\zeta \to 1/\zeta$.
Computing $\kappa(1/\zeta, w^2/X) \cdot \kappa(\zeta, X)$ symbolically, the parameter $\zeta$ cancels completely from the product:
$$\kappa(1/\zeta, w^2/X) \cdot \kappa(\zeta, X) = \frac{w^2 (s X^N - 1)(s (w^2/X)^N - 1)}{X^{N/2} (w^2/X)^{N/2} (w^2 - 1)^2} = \frac{s\cosh 2T\theta - \cos\pi T}{2\sin^2(\pi\omega/2)} .$$
Since every eigenvalue of $K(-\theta)^{\rm T} K(\theta)$ equals $U$, $K(-\theta)^{\rm T} K(\theta) \equiv U \mathbf 1$ identically for all $N$ and $\theta$. The proof is complete and rigorous. Confirmed as [Theorem].

### Claim 4
**Crossing function $\Phi_n$.**
The boundary crossing-unitarity relation $k(\theta) = \Phi \check R_{\bar V\bar V}(2\theta) k(-\theta)$ was evaluated for $n=1..5$ and $s=\pm 1$ across multiple values of $\omega \in \{0.37, 1.83, 2.6\}$.
The numerical check confirmed:
$$\Phi_+ = \frac{\sin\pi(\omega/2 - \mu)}{\sin\pi(\omega/2 + \mu)}, \qquad \Phi_- = \frac{\cos\pi(\mu - \omega/2)}{\cos\pi(\mu + \omega/2)}$$
with residual $\le 1.3 \times 10^{-14}$.
For $n=2$, with $T = 3\omega/2 = 3\lambda/2$, the resulting functional equations match Gandenberger (hep-th/9806003) eqs. (4.4) and (4.6). Confirmed.

### Claim 5
**Singlet projection.**
The operator $(1 \otimes K)\check R(1 \otimes K^{\rm T})$ acting on the $B_1$ singlet in $V \otimes \bar V$ was evaluated for $n=2, 3, 4, 5$. In both orderings, the projected vector is strictly proportional to the singlet, with residual $\le 1.7 \times 10^{-14}$.
The eigenvalues match the closed forms:
$$c_+ = -\frac{\cos\pi\mu \cos\pi(\mu + \omega/2)\cos\pi(\mu - T)}{\sin^2(\pi\omega/2)\cos\pi(\mu - \omega/2)}, \qquad c_- = \frac{\sin\pi\mu \sin\pi(\mu + \omega/2)\sin\pi(\mu - T)}{\sin^2(\pi\omega/2)\sin\pi(\mu - \omega/2)} .$$
Confirmed.

### Claim 6
**Scalar factors $A^{(\pm)}$.**
The proposed closed forms for $A^{(\pm)}$ satisfy the functional equations $A(\mu)A(-\mu) = 1/U$ and $A(T/2 - \mu)/A(T/2 + \mu) = F_{11}(2\mu)/\Phi(\mu)$ to $3.5 \times 10^{-8}$ across $n=1..5$.
- At $n=1$: $A^{(+)}$ reproduces the Ghoshal–Zamolodchikov / BPTT factor at the Neumann point, and $A^{(-)}$ reproduces GZ at $(\pi/2, i\pi\omega/2)$, with discrepancy $< 7.3 \times 10^{-8}$.
- At $n=2$: $A^{(+)}$ matches Gandenberger eq. (4.14) times $\sigma^+ = (\lambda/4)(5\lambda/4)$ (eq. 5.3) to $2.5 \times 10^{-9}$.
- For $A^{(-)}$: inspection of Gandenberger's TeX source `paprev.tex` reveals that Gandenberger defined $K^-$ with an explicit factor of $i$ in the diagonal entry $\kappa_-(\mu) = i \frac{\cos(\pi(\mu - \lambda/4))}{\sin(\pi\lambda/2)}$ (eq. kappa / line 720), whereas the book's $K(s=-1)$ matrix has real entries on the diagonal. The relative factor of $-i$ between the worker's $A^{(-)}$ and Gandenberger's $\mathcal{A}^-$ is therefore an exact consequence of this normalization choice. Confirmed.

### Claim 7
**Soliton bootstrap.**
The bootstrap condition for the fused antisoliton $V^{\otimes n} \to \bar V$ was verified relative to $K=\mathbf 1$. The ratio method cancels bulk Clebsch-Gordan embedding norms and bulk S-matrix factors identically. The bootstrap identity $Y(\mu) = P(\mu) \prod_{k=1}^n Y(\mu + T - k\omega)$ holds for $Y = \rho_+ \sigma$ ($s=+1$, $n=2..5$) and $Y = -\rho_-$ ($s=-1$, $n=2..4$) with residual $\le 4.7 \times 10^{-9}$ ($5.8 \times 10^{-9}$ at $n=5$). Confirmed.

### Claim 8
**CDD factor analysis.**
We independently tested the bootstrap condition with and without the CDD factor $\sigma_{\rm extra} = (T/2 - \omega/2)_I(T/2 + \omega/2)_I$:
- For $n=2$ (even), omitting $\sigma_{\rm extra}$ causes the bootstrap to fail completely with $O(1)$ error ($1.30$). Thus $\sigma_{\rm extra}$ is strictly forced by the bootstrap for even $n$.
- For $n=3$ (odd), the bootstrap holds with error $2.2 \times 10^{-10}$ both with and without $\sigma_{\rm extra}$. However, omitting $\sigma_{\rm extra}$ changes $K_1$ by $k_{h/2-1} = 1 + O(B)$ with a nonzero $O(B)$ coefficient, which contradicts Kim's one-loop $a_3^{(1)}$ Neumann calculation (`kim1995` eqs. 33, 39).
- The factor $\sigma_{\rm extra}$ introduces a pole at $\theta = i\pi(1/2 - 1/h)$ in the physical strip. Direct inspection of Gandenberger's $a_2^{(1)}$ solution confirms that this exact pole is present in his published amplitude (via $\sigma^+ = (\lambda/4)(5\lambda/4)$ at $\mu = \lambda/4$, $\theta = i\pi/6$), where its physical interpretation was also left open. Confirmed.

### Claim 9
**Main result: breather amplitude and continuation.**
The breather amplitude identity $A(\mu_a)F_{11}(2\mu)A(\mu_b)c_+(\mu) = (-n\omega/2)_I(n\omega/2 - 1/2)_I(1/2)_I$ was verified numerically for $n=1..5$ to high precision (relative error $6.6 \times 10^{-9}$ at $n=3$, $6.3 \times 10^{-8}$ at $n=4$).
Under the continuation $\omega \to -2/B$, the blocks continue to $(-n)(n+B/2)(-B/2)$.
We proved symbolically in SymPy that for all $n \ge 1$:
$$(-n)\left(n+\frac B2\right)\left(-\frac B2\right) \equiv \frac{(0)(h-1+\frac B2)}{(h-1)(\frac B2)} \equiv K^N_1 .$$
The resulting amplitudes $K^N_a$ satisfy all boundary scattering axioms (unitarity, crossing, bootstrap) and tend to 1 as $B \to 0$ for all $n \le 7$. Confirmed.

### Claim 10
**The $s=-1$ bonus.**
For $s=-1$, the breather amplitude evaluates to DG(3.11)$\times \sigma_1(\mu)$.
Under continuation, it yields $K^{(-)}_a = -(a)(h-a)K_a^{\rm unif}$.
For $n=2$, this reproduces Gandenberger eq. (6.4), $K^{(-)}_1 = -(1)(-1+B/2)(3-B/2)$.
The amplitudes satisfy all axioms for $n \le 7$ and tend to 1 as $B \to 0$. Confirmed.

### Claim 11
**Semiclassics.**
The semiclassical check $\frac{\beta^2}{4h}\frac{d}{d\theta}(-i\ln A^{(+)}) \to \frac{1}{2}\ln(1 - \sin^2(\pi/h)/\cosh^2\theta)$ was re-evaluated at $\theta = 0.6$ for $\omega \in \{40, 80, 160\}$. The relative errors decrease as $1/\omega$ (e.g. for $n=2$: 0.0214, 0.0108, 0.0055), confirming the classical Neumann time advance. Confirmed.

### Claim 12
**Proposed book text.**
The proposed text for `sec-neumann-duality`, `sec-an-nondiagonal-k`, and open problem 3 accurately reflects the mathematical findings, respects book conventions, correctly labels statements as [Numerical] for $n \le 5$, and documents the odd-$n$ proviso.

## Errors and overclaims not tied to a single claim
None. The worker explicitly separated proved statements (the unitarity function $U_n$, the block continuation identity) from numerical checks ($n \le 5$), transparently highlighted the odd-$n$ CDD proviso and the unexplained pole at $\theta = i\pi(1/2 - 1/h)$, and did not modify repository files.

## Prior art found
- **Gandenberger 1999b (hep-th/9911178):** Derived the non-diagonal $a_n^{(1)}$ reflection matrices, but explicitly left their scalar factors undetermined.
- **Regelskis–Vlaar (arXiv:1803.06491) & Appel–Vlaar (arXiv:2007.09218):** Classified solutions of the reflection equations from quantum symmetric pairs and $q$-Onsager algebras, but did not compute the scalar factor for the non-diagonal Toda K-matrices.
- **Literature search:** No published closed form for the non-diagonal $a_n^{(1)}$ scalar factors or the soliton-side derivation of the Neumann particle amplitudes for $n \ge 3$ exists in prior literature. The closed forms presented in this report represent new results.

## Proposed book text
The replacement texts proposed in Section 7 of the worker's report are approved for integration into `sec-neumann-duality` and `sec-an-nondiagonal-k`.

## Items for an Opus pass
1. **Physical interpretation of the strip pole (Claim 8):** Attack the physical-strip pole at $\theta = i\pi(1/2 - 1/h)$ in $A^{(+)}$. Determine whether this pole can be interpreted as a boundary bound state or explained via a boundary Coleman–Thun diagram, following up on Gandenberger's open question for $n=2$.
2. **Odd-$n$ CDD factor mechanism (Claim 8):** Investigate whether a higher-representation bootstrap (e.g. bootstrapping from tensor products beyond the fundamental representation) or a higher-order boundary charge constraint uniquely forces the factor $(\frac T2 - \frac\omega2)_I(\frac T2 + \frac\omega2)_I$ for odd $n \ge 5$ without invoking uniformity in $n$.
3. **Analytic proof of the crossing function (Claim 4):** Attempt an analytic derivation of the crossing eigenvalue $\Phi_+$ for general $n$ from the spectral decomposition of $\check R_{\bar V\bar V}(2\theta)$.

## Scripts
- `/tmp/halfline-referee/03/check_all.py`: Master test suite checking claims 1–9.
- `/tmp/halfline-referee/03/particle_checks.py`: Verification of particle axioms, continuation, and $B \to 0$ limits.
- `/tmp/halfline-referee/03/paprev.tex`: Gandenberger's arXiv TeX source for exact normalization comparison.
- `/tmp/halfline-referee/03/test_cdd_bootstrap.py`: Independent check of the soliton bootstrap failure when omitting the CDD factor for even $n$.

Referee pass complete.
