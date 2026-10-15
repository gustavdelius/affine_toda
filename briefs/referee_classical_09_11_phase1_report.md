# Phase-1 referee report: Part III (chapters 9–11)

- **Date:** 2026-10-15
- **Model:** Gemini 4 Argon (High)
- **Commit reviewed:** `a03edd4896c927253721352a5bb253521c7f0f39`

---

## 1. Ranked risk list (written before checking)

1. **Rank 1 — `part3-classical/11-reality.qmd:32–41` (`#prp-hirota-energy`) and `11-reality.qmd:45–47`: Unconditional statement of the energy–momentum formula $E=\sum_i M_{a_i}\cosh\theta_i$, $P=\sum_i M_{a_i}\sinh\theta_i$ for constituent solitons with complex rapidities, and its citations to `[@tw1999, eq. 4.6]` and `[@khastgir1996]`.**
   - *Claim quoted:* `"For a Hirota solution, the energy and momentum densities are total derivatives [@olive1993a; @olive1993b]. They integrate to $E=\sum_iM_{a_i}\cosh\theta_i,\qquad P=\sum_iM_{a_i}\sinh\theta_i,$ where $(a_i,\theta_i)$ are the species and (possibly complex) rapidities of the constituent solitons [@tw1999, eq. 4.6]."` and `"The complex-energy solutions of Khastgir and Sasaki have asymptotic data that are not $\Theta$-invariant."`
   - *Why it is at risk:* Evaluating $[\partial_x\ln\tau_j]_{-\infty}^{+\infty}$ produces $\sum_i m_{a_i}\cosh\theta_i$ only when $\operatorname{Re}(\cosh\theta_i)>0$ for every constituent (if $\operatorname{Re}(\cosh\theta_i)<0$, the roles of $x\to+\infty$ and $x\to-\infty$ swap for $E_i$, flipping the sign of its contribution unless continued analytically from real rapidities), when the multi-soliton interaction coefficient $\prod_{k<l}A_{a_ka_l}(\theta_k-\theta_l)$ is neither $0$ nor $\infty$, and when no $\tau_j$ vanishes on the integration contour. Furthermore, Khastgir and Sasaki (`[@khastgir1996]`) also discuss negative-energy solutions (`10-solitons.qmd:180`), whose relation to $\Theta$-invariance and `@prp-hirota-energy` is not explained.
   - *How I will check it:* Read `hep-th/9810006` (`tw1999`), `hep-th/9507001` (`khastgir1996`), and `hep-th/9305160` (`olive1993b`), and write an independent script in `referee_scratch/classical_09_11/` computing the integrated energy and momentum of 1-soliton and 2-soliton Hirota solutions at complex rapidities (including $\operatorname{Re}\cosh\theta_i<0$ and zeros of $A_{ab}$).

2. **Rank 2 — `part3-classical/11-reality.qmd:43–49` (`#sec-reality-energy-momentum`) and `appendices/E-original-results.qmd:8`: The $\Theta$-selection rule (the sole claimed original result of Part III), its hidden hypotheses, its status label in Appendix E versus Chapter 11, and the $a_2^{(1)}$ spin-2 example.**
   - *Claim quoted:* `"If the asymptotic data are $\Theta$-invariant, i.e. the multiset $\{(a_i,\theta_i)\}$ is mapped to itself by $(a,\theta)\mapsto(\bar a,\bar\theta)$, then every conserved charge satisfies $Q=\epsilon_Q\bar Q$, where $J_Q(-\bar{\vec\phi})=\epsilon_Q\overline{J_Q(\vec\phi)}$ defines $\epsilon_Q=\pm1$; in particular $E$ and $P$ are real."`, `11-reality.qmd:46`: `"In $a_2^{(1)}$ the spin-2 charge of a single soliton is $q_ae^{2\theta}$, where by Freeman's computation [@freeman1995] the $q_a$ form the eigenvector of the $A_2$ Cartan matrix with eigenvalue $3$, so $q_2=-q_1$, and they can be normalized to be real. If the charge is additive over the constituents, as $E$ and $P$ are, a breather $\{(1,\theta_0+iu),(2,\theta_0-iu)\}$ has spin-2 charge $q_1\big(e^{2(\theta_0+iu)}-e^{2(\theta_0-iu)}\big)=2iq_1e^{2\theta_0}\sin2u$, which is imaginary. This agrees with the selection rule with $\epsilon_Q=-1$."`, and `appendices/E-original-results.qmd:8`: `"- **[Theorem, given additivity of the higher charges]** A selection rule from the antilinear symmetry $\Theta$... (@sec-reality-energy-momentum)."`
   - *Why it is at risk:* (a) In `11-reality.qmd`, the selection rule is an unlabelled bullet point, whereas Appendix E gives it the status label `**[Theorem, given additivity of the higher charges]**`. (b) For a solution whose $(a_i,\theta_i)$ multiset is $\Theta$-invariant (such as a breather with $\xi_2=\bar\xi_1+i\pi$), the field $\vec\phi$ itself is not $\Theta$-invariant ($\Theta\vec\phi\ne\vec\phi$), so concluding $Q[\vec\phi]=Q[\Theta\vec\phi]=\epsilon_Q\overline{Q[\vec\phi]}$ requires not only additivity of $Q$, but also independence of $Q$ from the position/phase parameters $\xi_i$, absence of contour-crossing singularities between $\xi_i$ and $\bar\xi_i$, and that $\Theta:\tau_j\mapsto\bar\tau_j$ acts on constituent data by $(a,\theta,\xi)\mapsto(\bar a,\bar\theta,\bar\xi)$ in general algebras (not only $a_n^{(1)}$). (c) $J_Q(-\bar{\vec\phi})=\epsilon_Q\overline{J_Q(\vec\phi)}$ holds only in a $\Theta$-eigenbasis of charges (relevant when an exponent has multiplicity $>1$, as in $d_{2n}^{(1)}$) and depends on the overall phase normalization of $J_Q$; we must check Freeman (`hep-th/9408092`) and explicitly construct the $a_2^{(1)}$ spin-2 density to verify $\epsilon_Q$ and $q_a$.
   - *How I will check it:* Read `hep-th/9408092` (`freeman1995`), construct the $a_2^{(1)}$ spin-2 conserved density explicitly in an independent script, evaluate it on single solitons and breathers, and check how $\Theta:\tau_j\mapsto\bar\tau_j$ acts on single solitons of $d_n^{(1)}$ and $e_6^{(1)}$.

3. **Rank 3 — `part3-classical/10-solitons.qmd:128–136` (`#sec-breathers`), `11-reality.qmd:28` (`#fig-breather-theta`), and `11-reality.qmd:87` (`#exr-breather-real`): Parameter range $u\in(0,\pi/2)$, interaction coefficient $A_{a,\bar a}(2iu)$, and singularity conditions for $a_n^{(1)}$ breathers.**
   - *Claim quoted:* `"A soliton of species $a$ and one of the conjugate species $\bar a=h-a$, with complex-conjugate rapidities $\theta_0\pm iu$ (so $u$ here is $\frac\pi2$ minus the $u$ of @sec-sg-classical) and suitably related $\xi$'s, form a breather: a localized, time-periodic solution of zero topological charge."`, `11-reality.qmd:28`: `"Right: the classical breather mass $2M_a\cos u$, for $u\in(0,\pi/2)$, fills a continuous band."`, and `11-reality.qmd:87`: `"If moreover $\xi_2=\bar\xi_1$, show that $\bar\tau_j=\tau_j$, so that the field is purely imaginary where all $\tau_j>0$, and that the $\tau_j$ change sign, so this solution is singular along curves in the $(x,t)$ plane."`
   - *Why it is at risk:* For $a_n^{(1)}$ with $h\ge 3$ and $a\ne h/2$, $A_{a,h-a}(2iu)=\frac{\cos 2u+\cos(2\pi a/h)}{\cos 2u+1}$ vanishes at $u=\frac{\pi}{2}-\frac{\pi a}{h}\in(0,\pi/2)$ and is positive for $0<u<\frac{\pi}{2}-\frac{\pi a}{h}$. At $u=\frac{\pi}{2}-\frac{\pi a}{h}$ the $E_1E_2$ term in $\tau_j$ disappears, so the solution is not a localized breather of zero topological charge and mass $2M_a\cos u$. Furthermore, in Harder, Iskandar and McGhee (`hep-th/9409035`), we must check the allowed range of $u$ for a genuine breather and check whether $\tau_j$ changes sign for $\xi_2=\bar\xi_1$ across $(0,\pi/2)$.
   - *How I will check it:* Read `hep-th/9409035` (`harder1995`) and analyze $\tau_j(x,t)$ for $a_n^{(1)}$ breathers analytically and numerically in `referee_scratch/classical_09_11/`.

4. **Rank 4 — `part3-classical/10-solitons.qmd:174` (`#sec-classical-literature`): Pointwise total-derivative formula for the energy–momentum tensor $T_{\mu\nu}=(\eta_{\mu\nu}\partial^2-\partial_\mu\partial_\nu)C$.**
   - *Claim quoted:* `"Although the fields are complex, on Hirota solutions the energy–momentum tensor is a total derivative, $T_{\mu\nu}=(\eta_{\mu\nu}\partial^2-\partial_\mu\partial_\nu)C$ with $C=-\frac2{\beta^2}\sum_j\frac2{\vec\alpha_j^2}\ln\tau_j$, and the energies are real [@olive1993a; @olive1993b; @zhu1995]."`
   - *Why it is at risk:* While the trace $T^\mu{}_\mu=2V_{\rm imag}=\partial^2 C$ follows directly from the Hirota bilinear equations (`eq-hirota-bilinear`), the individual components $T_{00},T_{01},T_{11}$ of the canonical stress tensor $T_{\mu\nu}=\partial_\mu\vec\phi\cdot\partial_\nu\vec\phi-\eta_{\mu\nu}\mathcal L$ need not equal $(\eta_{\mu\nu}\partial^2-\partial_\mu\partial_\nu)C$ pointwise unless either an improvement term is added to $T_{\mu\nu}$ or a non-trivial bilinear identity holds for all Hirota solutions.
   - *How I will check it:* Read `hep-th/9305160` (`olive1993b`) and `hep-th/9404097` (`zhu1995`), and test pointwise whether $T_{\mu\nu}^{\rm can}=(\eta_{\mu\nu}\partial^2-\partial_\mu\partial_\nu)C$ on 1-soliton and 2-soliton solutions of $a_1^{(1)}, a_2^{(1)}, a_3^{(1)}$ in an independent script.

5. **Rank 5 — `part3-classical/10-solitons.qmd:22–25` (`#sec-imaginary-coupling-vacua`), `10-solitons.qmd:104` (`#sec-soliton-masses`), and `10-solitons.qmd:158–160` (`#exr-lie-duality`): Soliton masses for twisted and untwisted non-simply-laced affine algebras and Lie duality.**
   - *Claim quoted:* `"The soliton masses are proportional to particle masses: - in the simply-laced and twisted cases, to the particle masses of the same theory; - in the untwisted non-simply-laced cases, to those of the theory based on $(g^\vee)^{(1)}$ ("Lie duality")."` and `10-solitons.qmd:104`: `"For the other untwisted non-simply-laced algebras it is quoted from [@kneipp1996]."`
   - *Why it is at risk:* In folded theories, a single soliton of the folded theory is obtained by folding a parent soliton or an orbit of parent solitons, and the mass ratios depend on how the parent species fold (and on whether the parent solitons in an orbit interact or not). We must check Kneipp and Olive (`hep-th/9404030`) and McGhee (`hep-th/9307035`) to see whether the soliton masses in all twisted cases ($a_{2n-1}^{(2)}, d_{n+1}^{(2)}, e_6^{(2)}, d_4^{(3)}, a_{2n}^{(2)}$) and untwisted non-simply-laced cases ($b_n^{(1)}, c_n^{(1)}, f_4^{(1)}, g_2^{(1)}$) obey the stated rule, and what the proportionality constant $M_a/m_a$ is.
   - *How I will check it:* Read `hep-th/9404030` (`kneipp1996`) and `hep-th/9307035` (`mcghee1994`), and compute the folded soliton masses and particle masses for all folded families in an independent script.

6. **Rank 6 — `part3-classical/10-solitons.qmd:76–77` (`#sec-vertex-solitons`): Polynomial degrees of $\tau_j$ in folded theories and the explicit $c_2^{(1)}$ folded soliton $\tau$-functions.**
   - *Claim quoted:* `"For folded theories the relevant label is $n_j^\vee$, the level of $L(\Lambda_j)$. A coincident orbit $O$ of $\lvert O\rvert$ parent solitons has degree $\lvert O\rvert n_j^\vee$. In $c_n^{(1)}$ every $n_j^\vee$ is 1, so a coincident pair of species $a$ and $2n-a$ has degree 2 at every node, including the long end nodes, whose Kac label is 1, and the species-$n$ soliton has degree 1 at every node, including the short nodes, whose Kac label is 2. For $c_2^{(1)}$, folded from $a_3^{(1)}$, the static pair of species 1 and 3 has $\tau_0=1+2E+\frac12E^2$, $\tau_1=\tau_3=1+\frac12E^2$ and $\tau_2=1-2E+\frac12E^2$."`
   - *Why it is at risk:* (a) In $c_2^{(1)}$ there are only three nodes ($0,1,2$), so writing $\tau_1=\tau_3$ uses the parent $a_3^{(1)}$ indexing without stating the $c_2^{(1)}$ indexing. (b) We must check whether $\tau_0=1+2E+\frac12E^2$, $\tau_1=1+\frac12E^2$, $\tau_2=1-2E+\frac12E^2$ actually satisfies the $c_2^{(1)}$ bilinear equations (`eq-hirota-bilinear`) with $n_0=1, n_1=2, n_2=1$ and $\alpha_0^2=\alpha_2^2=2, \alpha_1^2=1$. (c) In $c_n^{(1)}$ the parent theory is $a_{2n-1}^{(1)}$, where every parent representation has level $n_j^{\rm parent}=1$; in general foldings (e.g. $d_{n+1}^{(1)}\to b_n^{(1)}$ or $d_4^{(1)}\to g_2^{(1)}$), we must check whether the degree of $\tau_j$ for a folded orbit $O$ of parent solitons is $|O|n_j^{\rm parent}$ rather than $|O|n_j^\vee$.
   - *How I will check it:* Check the $c_2^{(1)}$ bilinear equation symbolically, and compare the degree formula $|O|n_j^\vee$ vs $|O|n_j^{\rm parent}$ across all foldings in `@tbl-foldings` (`b_n^{(1)}, c_n^{(1)}, f_4^{(1)}, g_2^{(1)}, a_{2n-1}^{(2)}, d_{n+1}^{(2)}, e_6^{(2)}, d_4^{(3)}$).

7. **Rank 7 — `part3-classical/11-reality.qmd:58–62` (`#sec-singular-solutions`) and `11-reality.qmd:90–92` (`#exr-singular-points`): Zeros of $\tau_j$ in the real $(x,t)$ plane, double poles with zero residue, and whether integrated charges are unaffected.**
   - *Claim quoted:* `"Zeros of $\tau_j$ in the real $(x,t)$ plane are generically isolated points: one complex condition, two real conditions. Non-generic data, such as a single soliton with tuned $\operatorname{Im}\xi$ or a breather with $\xi_2=\bar\xi_1$, give curves of zeros. There the energy density has double poles with zero residue, so integrated charges are unaffected, and the Hirota data continue analytically through the singularity."` and `@exr-singular-points`: `"Show that a zero of $\tau_j$ in the real $(x,t)$ plane is generically an isolated point, and that near it the energy density has a double pole with zero residue, so that integrated charges are unaffected."`
   - *Why it is at risk:* (a) For a single soliton with real rapidity $\theta$, $\tau_j(x,t)$ is a function of a single real linear combination $x\cosh\theta-t\sinh\theta$, so its real zero set is never an isolated point in $\mathbb R^2$ (it is empty for generic $\operatorname{Im}\xi$, and a line for tuned $\operatorname{Im}\xi$). (b) When a zero of $\tau_j$ lies on the real $x$-axis at time $t$, a double pole $c/(x-x_0)^2$ is not Lebesgue/improper-Riemann integrable on $\mathbb R$: it integrates to the finite value only via contour deformation into $\mathbb C$ (which is independent of whether the contour passes above or below $x_0$ because the residue vanishes). (c) The topological charge $\frac{\beta}{2\pi}\big(\vec\phi(+\infty)-\vec\phi(-\infty)\big)$ is *not* unaffected when a zero of $\tau_j$ crosses the real $x$-axis: it jumps by a coweight/coroot!
   - *How I will check it:* Compute the Laurent expansion of $T_{00}^{\rm can}$ and higher charge densities around a simple zero $x_0(t)$ of $\tau_j(x,t)$, and check the behaviour of the energy and topological charge across a zero crossing.

8. **Rank 8 — `part3-classical/09-lagrangian.qmd:90–96` (`#eq-spin3-current`) and `09-lagrangian.qmd:168–170` (`#exr-spin3`): The spin-3 conserved current of sinh-Gordon and sine-Gordon, and its evaluation on the sine-Gordon kink.**
   - *Claim quoted:* `"\partial_-\Big[(\partial_+^2\phi)^2+\frac{\beta^2}4(\partial_+\phi)^4\Big]+\partial_+\Big[m^2\cosh(\beta\phi)\,(\partial_+\phi)^2\Big]=0 . ... The charge $Q_3=\int dx\,\big(T_4+\mathcal X_2\big)$ built from this current has spin 3, and on a one-particle state it takes the value $\chi\,m^3e^{3\theta}$, as in (-@eq-charge-eigenvalue). For sine-Gordon, replace $\beta^2$ by $-\beta^2$."`
   - *Why it is at risk:* Passing from sinh-Gordon ($\partial_+\partial_-\phi=-\frac{m^2}{\beta}\sinh\beta\phi$) to sine-Gordon ($\partial_+\partial_-\phi=-\frac{m^2}{\beta}\sin\beta\phi$) requires $\beta\to i\beta$, which replaces $\frac{\beta^2}{4}$ by $-\frac{\beta^2}{4}$ *and* replaces $\cosh(\beta\phi)$ by $\cos(\beta\phi)$ (and subtracts 1 if integrating on a kink when $\phi(\pm\infty)$ is a non-zero vacuum, though $\partial_+\phi\to 0$ at $\pm\infty$). Literally replacing only $\beta^2\to -\beta^2$ leaves $\cosh(\beta\phi)$, which violates conservation for sine-Gordon!
   - *How I will check it:* Differentiate the current symbolically using the sine-Gordon and sinh-Gordon equations of motion, and evaluate $Q_3$ on the sine-Gordon kink (`eq-sg-kink`).

9. **Rank 9 — `part3-classical/10-solitons.qmd:82–90` (`#sec-topological-charges`), `10-solitons.qmd:92` (`#fig-topological-charges`), and `10-solitons.qmd:154–156` (`#exr-charges-a3`): Topological charges of $a_n^{(1)}$ single solitons as weights of $\Lambda^a\mathbb C^h$, and the jumping rule as $\operatorname{Im}\xi$ varies.**
   - *Claim quoted:* `"\vec\phi(+\infty)=-\frac1\beta\sum_j\vec\alpha_j\,\big[\vartheta_j\big]_{(-\pi,\pi)} . This lies in $\frac{2\pi}\beta\Lambda_W^\vee$. The topological charge $\frac\beta{2\pi}\vec\phi(+\infty)$ is a weight of the $a$-th fundamental representation $\Lambda^a\mathbb C^h$ of $su(h)$ [@mcghee1994]. Which weight depends on $\operatorname{Im}\xi$: as $\operatorname{Im}\xi$ increases through $2\pi$, the reduced phases jump at $h/\gcd(a,h)$ values, $\gcd(a,h)$ of them at once (those with $j$ in one residue class mod $h/\gcd(a,h)$), and the charge runs through a cycle of $h/\gcd(a,h)$ weights."`
   - *Why it is at risk:* Sign conventions in `eq-hirota-ansatz` ($\vec\phi=+\frac{i}{\beta}\sum_j\vec\alpha_j^\vee\ln\tau_j$) and `eq-an-soliton` ($\omega=e^{+2\pi i/h}$) determine whether $\frac{\beta}{2\pi}\vec\phi(+\infty)=-\frac{1}{2\pi}\sum_{j=0}^{h-1}\vec\alpha_j[\vartheta_j]_{(-\pi,\pi)}$ lies in the fundamental representation $\Lambda^a\mathbb C^h$ (highest weight $\vec\lambda_a$) or the conjugate representation $\Lambda^{h-a}\mathbb C^h$ (highest weight $\vec\lambda_{h-a}$).
   - *How I will check it:* Compute $-\frac{1}{2\pi}\sum_{j=0}^{h-1}\vec\alpha_j[\vartheta_j]_{(-\pi,\pi)}$ in terms of the fundamental weights $\vec\lambda_1,\dots,\vec\lambda_n$ for $a_2^{(1)}, a_3^{(1)}, a_4^{(1)}$ and all $a=1,\dots,n$ in an independent script, and check `hep-th/9307035` (`mcghee1994`).

10. **Rank 10 — `part3-classical/09-lagrangian.qmd:84–88` (`#sec-conserved-charges`): Spins of conserved charges (exponents modulo $kh$) for all untwisted and twisted affine Kac–Moody algebras, and the claim that dual pairs have the same spins.**
    - *Claim quoted:* `"s\in\{\text{exponents of }\hat g\}+kh\,\mathbb Z , where $h=\sum_jn_j$ and $kh$ is the order of the twisted Coxeter element ([@kac1990] lists the exponents of $X_N^{(k)}$ modulo $kh$). For untwisted $\hat g=g^{(1)}$ these are the exponents of $g$ (@tbl-simple-lie) repeated modulo $h$. For $a_2^{(2)}$ they are $1,5$ mod 6; for $d_4^{(3)}$, $1,5,7,11$ mod 12; for $e_6^{(2)}$, $1,5,7,11,13,17$ mod 18; for $a_{2n}^{(2)}$, the odd integers not divisible by $2n+1$; for $a_{2n-1}^{(2)}$ and $d_{n+1}^{(2)}$ (as for $b_n^{(1)}$, $c_n^{(1)}$), all odd integers. The two members of a dual pair have the same spins."`
    - *Why it is at risk:* With the book's definition $h=\sum_j n_j$ (sum of coprime Kac labels), we must verify whether $kh$ is indeed the order of the twisted Coxeter element for every twisted algebra $X_N^{(k)}$ ($a_{2n}^{(2)}, a_{2n-1}^{(2)}, d_{n+1}^{(2)}, e_6^{(2)}, d_4^{(3)}$) and whether the explicit lists of exponents and their agreement across dual pairs are completely accurate.
    - *How I will check it:* Compute the principal degrees of the centralizer of $E_1=\sum_j c_j e_j$ for all untwisted and twisted affine algebras directly in an independent script.

11. **Rank 11 — `part3-classical/09-lagrangian.qmd:122–138` (`#sec-masses-couplings`), `#eq-area-rule`, and `#exr-shg-cubic` (`09-lagrangian.qmd:180–182`): Three-point couplings $C_{abc}$, specification of the charge eigenbasis for degenerate masses ($a_n, d_n, e_6$), and the area rule $|C_{abc}|=\frac{4\beta}{\sqrt{h}}\Delta_{abc}$.**
    - *Claim quoted:* `"For degenerate masses ($a_n$, $d_n$, $e_6$) the $\vec e_a$ must be charge eigenvectors, the complex eigenvectors of the centre of the affine diagram (of $\mathbb Z_h$ for $a_n$), normalized by $\vec e_a\cdot\vec e_{\bar a}=1$ so that $\phi_{\bar a}=\bar\phi_a$. Particle $a$ is the eigenvector with the centre character of node $a$. In a real orthonormal basis of a degenerate pair neither (-@eq-area-rule) nor @prp-dorey-rule holds. For simply-laced $\hat g$ these have a striking geometric form [@braden1990]: $\lvert C_{abc}\rvert=\frac{4\beta}{\sqrt h}\,\Delta_{abc}$."`
    - *Why it is at risk:* (a) "the centre of the affine diagram" and "the centre character of node $a$" are imprecise for $d_n^{(1)}$ and $e_6^{(1)}$, where the degenerate particles are conjugate pairs under the *diagram involution of the finite Dynkin diagram* $\sigma$ (which fixes the affine node $0$), not under the centre of $G$ (which moves node $0$), and for even $n$ the two spinor particles of $d_n^{(1)}$ are each self-conjugate ($\bar s=s, \bar s'=s'$) and are distinguished by the $\mathbb Z_2$ diagram flip of the spinor fork! (b) We must also check whether in $d_{2m}^{(1)}$ a real orthonormal basis of $(s,s')$ — namely the eigenvectors of the fork flip — *is* the physical basis in which `eq-area-rule` and `@prp-dorey-rule` hold!
    - *How I will check it:* Compute $\mathsf M^2$, the diagram automorphisms, $C_{abc}$, and $\frac{4\beta}{\sqrt{h}}\Delta_{abc}$ for $a_2^{(1)}, a_3^{(1)}, d_4^{(1)}, d_5^{(1)}, e_6^{(1)}, e_7^{(1)}, e_8^{(1)}$ in an independent script.

12. **Rank 12 — `part3-classical/09-lagrangian.qmd:142–160` (`#sec-dorey-rule`, `#prp-dorey-rule`, `#fig-dorey-d4`, `#exr-dorey-a2`): Bicoloured Coxeter element $w=s_{\rm black}s_{\rm white}$, representative roots $\vec\gamma_a=(1-w^{-1})\vec\lambda_a$, and the $D_4$ triangle in `@fig-dorey-d4`.**
    - *Claim quoted:* `"If $a$ is white, the black reflections fix $\vec\lambda_a$ and $s_{\rm white}\vec\lambda_a=\vec\lambda_a-\vec\alpha_a$, so $\vec\gamma_a=\vec\alpha_a$. If $a$ is black, the same argument gives $\vec\gamma_a=s_{\rm white}\vec\alpha_a$."` and `@fig-dorey-d4`: `"The three roots $\vec\gamma_1+\vec\gamma_1'+\vec\gamma_c=0$ close a triangle with sides $1,1,\sqrt3$ and angles $\pi/6,\pi/6,2\pi/3$: two outer particles fuse into the central one."`
    - *Why it is at risk:* In $d_4^{(1)}$, we must verify whether two roots from the *same* outer orbit $\Omega_1$ sum to a root in the central orbit $-\Omega_c$ (so $C_{11c}\ne 0$), and check the SVG figure `part3-classical/figures/ch9_f13_dorey_d4.svg` and the derivations of $\vec\gamma_a$ and `@exr-dorey-a2`.
    - *How I will check it:* Compute the orbits $\Omega_a$ and all root triangles for $A_2$ and $D_4$ in an independent script and inspect `ch9_f13_dorey_d4.svg`.

13. **Rank 13 — `part3-classical/10-solitons.qmd:122` (`#sec-multi-solitons`): Classical spatial shift $\Delta x$, time delay $\Delta t$, Wigner's relation $\Delta t=d\delta/dE_1$ in the rest frame of soliton 2, and the leading semiclassical phase of sine-Gordon $S_0(\theta)$.**
    - *Claim quoted:* `"Since $M_{a_1}=2hm_{a_1}/\beta^2$ (-@eq-soliton-mass), the shift is $E_1\Delta x=-p_1\Delta t=-(2h/\beta^2)\ln A_{a_1a_2}$, with $E_1,p_1$ the energy and momentum of soliton 1 [@fring1994]. Semiclassically, write the S-matrix as $e^{i\delta_{ab}(\theta)}$. Wigner's relation $\Delta t=d\delta/dE_1$, in the rest frame of soliton 2, gives $d\delta_{ab}/d\theta=p_1\Delta t=(2h/\beta^2)\ln A_{ab}(\theta)$ at leading order. It is negative, as for an advance. For sine-Gordon this is the leading term of $-i\,d\ln S_0/d\theta$ from (-@eq-sg-s0), namely $(16/\beta_{\rm SG}^2)\ln\tanh\frac\theta2$."`
    - *Why it is at risk:* Sign and frame conventions for $\Delta x, \Delta t, d\delta/d\theta$, the small-$\xi$ limit of `eq-sg-s0`, and the reuse of $E_1$ for both the exponential of soliton 1 and the energy of soliton 1 in the same paragraph.
    - *How I will check it:* Read `hep-th/9405034` (`fring1994`), derive the asymptotic worldlines as $t\to\pm\infty$, and compute the small-$\xi$ asymptotics of `-i d\ln S_0/d\theta` from `eq-sg-s0`.

14. **Rank 14 — `part3-classical/09-lagrangian.qmd:114–117` (`#sec-masses-couplings`) and `#exr-dn-masses` (`09-lagrangian.qmd:172–174`): Explicit mass formulas for $a_n^{(1)}$, $d_n^{(1)}$ and $e_8^{(1)}$, and the sum rule $\sum_a m_a^2=m^2\sum_j n_j\vec\alpha_j^2$.**
    - *Claim quoted:* Table in `09-lagrangian.qmd:114–117` giving $m_a$ for $a_n^{(1)}$, $d_n^{(1)}$, and trigonometric expressions for $m_1,\dots,m_8$ of $e_8^{(1)}$ in terms of $m_1$ and $m_2$, plus $\sum_a m_a^2=m^2\sum_j n_j\vec\alpha_j^2=2hm^2$.
    - *Why it is at risk:* The eight closed-form trigonometric expressions for $e_8^{(1)}$ and the normalization of $d_n^{(1)}$ ($2\sqrt2 m\sin\frac{\pi a}{h}$ and $\sqrt2 m$) must match the eigenvalues of $\mathsf M^2=m^2\sum_j n_j\vec\alpha_j\vec\alpha_j^{\mathsf T}$ exactly.
    - *How I will check it:* Compute the eigenvalues of $\mathsf M^2$ symbolically/numerically to 50 digits in an independent script and compare with every formula in the table.

15. **Rank 15 — `part3-classical/10-solitons.qmd:10`: Embedded sine-Gordon kinks and the species-$n$ soliton of $a_{2n-1}^{(1)}$ at $\operatorname{Im}\xi=\pi/2$.**
    - *Claim quoted:* `"(Exceptions are embedded sine-Gordon kinks, real solutions of the form $\vec\phi=\vec v\,\varphi(x,t)$ with a fixed vector $\vec v$ and $\varphi$ a sine-Gordon kink. An example is the species-$n$ soliton of $a_{2n-1}^{(1)}$ at $\operatorname{Im}\xi=\pi/2$, which has $\tau_{2k+1}=\bar\tau_{2k}$.)"`
    - *Why it is at risk:* Must check that $\vec\phi(x,t)$ is strictly real, satisfies $\tau_{2k+1}=\bar\tau_{2k}$, and reduces to a sine-Gordon kink along a fixed direction $\vec v\in\mathbb R^{2n-1}$.
    - *How I will check it:* Derive $\tau_j$ and $\vec\phi$ analytically for species $n$ of $a_{2n-1}^{(1)}$ at $\operatorname{Im}\xi=\pi/2$ and verify in a script.

16. **Rank 16 — `part3-classical/09-lagrangian.qmd:57–79` (`#sec-lax-pair`, `#eq-lax-pair`, `#eq-zero-curvature`, `#exr-lax`): The Lax pair, zero-curvature equation, Lorentz boost transformation of $\mathcal A_\pm$, and `@exr-lax`.**
    - *Claim quoted:* `eq-lax-pair`, `eq-zero-curvature`, line 77: `"\mathcal A_\pm[\vec\phi_\vartheta](x;\lambda)=e^{\mp\vartheta}\mathcal A_\pm[\vec\phi](x';e^{\vartheta}\lambda)"`, and `@exr-lax`.
    - *Why it is at risk:* Sign and normalization conventions in $\partial_\pm=\partial_t\pm\partial_x$, $c_j=\sqrt{n_j\vec\alpha_j^2/2}$, the boost transformation $\vec\phi_\vartheta(x^+,x^-)=\vec\phi(e^{-\vartheta}x^+,e^{\vartheta}x^-)$, and the $2\times 2$ matrix representation for $a_1^{(1)}$ in `@exr-lax`.
    - *How I will check it:* Re-derive `eq-zero-curvature`, the boost identity, and `@exr-lax` by hand and symbolically in Python.

17. **Rank 17 — `part3-classical/10-solitons.qmd:68–79` (`#sec-vertex-solitons`) vs `part2-algebra/06-affine-algebras.qmd:165–171`: The vertex-operator formula for $\tau_j$ and the nilpotency bound $(F^a)^{n_j+1}=0$.**
    - *Claim quoted:* `"\tau_j=\langle\Lambda_j\rvert\,\prod_k\exp\big(Q_kF^{a_k}(z_k)\big)\,\lvert\Lambda_j\rangle ... For simply-laced $\hat g$, $(F^a)^{n_j+1}=0$ on $L(\Lambda_j)$, which has level $n_j$ [@olive1993b], so a single-soliton $\tau_j$ has degree at most $n_j$ in $E$."`
    - *Why it is at risk:* Compare the omission of $W_k(x,t)$ in `10-solitons.qmd:70` with `06-affine-algebras.qmd:168`, and verify the level and nilpotency proof in Olive, Turok and Underwood (`hep-th/9305160`).
    - *How I will check it:* Read `hep-th/9305160` (`olive1993b`) and compare with `10-solitons.qmd:68–79`.

18. **Rank 18 — `part3-classical/10-solitons.qmd:27` (`#sec-imaginary-coupling-vacua`): The quantum group $U_q(\hat g^\vee)$ and its deformation parameter $q=e^{4\pi^2 i/\beta^2}$ across all affine algebras.**
    - *Claim quoted:* `"The quantum theory has a non-local $U_q(\hat g^\vee)$ symmetry with $q=e^{4\pi^2i/\beta^2}$ for the coproduct of @sec-quantum-groups, and $q^{-1}$ for the opposite one used in Part VI (-@eq-q-physical) [@bernard1991; @felder1992]."`
    - *Why it is at risk:* Compare the sign and definition of $q$ in `10-solitons.qmd:27` with `eq-q-physical` in Chapter 20 (`part6-soliton-smatrices/20-quantum-group-symmetry.qmd`), where $q=e^{-4\pi^2 i/\beta^2}=-e^{-i\pi\omega}$ is fixed for Part VI's coproduct, and check whether non-simply-laced or twisted algebras modify the exponent.
    - *How I will check it:* Inspect `part2-algebra/07-quantum-groups.qmd` and `part6-soliton-smatrices/20-quantum-group-symmetry.qmd` and check Bernard–LeClair (`@bernard1991`).

19. **Rank 19 — `part3-classical/10-solitons.qmd:126` (`#sec-multi-solitons`) and `10-solitons.qmd:162–164` (`#exr-fusion`): Classical soliton fusion in $a_n^{(1)}$ at $\theta_1-\theta_2=i\pi(a+b)/h$.**
    - *Claim quoted:* `"At $\theta_1-\theta_2=i\pi(a+b)/h$, the denominator of $A_{ab}$ vanishes. Take $\theta_1=\theta+i\pi b/h$ and $\theta_2=\theta-i\pi a/h$. Then $E_1E_2$ has the exponent of a single soliton of species $a+b$ with rapidity $\theta$, because $m_ae^{i\pi b/h}+m_be^{-i\pi a/h}=m_{a+b}$, and the phases combine to $\omega^{j(a+b)}$. Letting $\xi_1,\xi_2\to-\infty$ while keeping $A_{ab}e^{\xi_1+\xi_2}$ fixed reduces (-@eq-n-soliton) to the single soliton of species $a+b$."`
    - *Why it is at risk:* Setting $\theta_1-\theta_2=i\pi(a+b)/h$ before taking $\xi_1,\xi_2\to-\infty$ makes $A_{ab}=\infty$, whereas the limit requires $\theta_1-\theta_2\to i\pi(a+b)/h$ with $\xi_1,\xi_2\to-\infty$ jointly so that $\ln A_{ab}(\theta_1-\theta_2)+\xi_1+\xi_2\to\xi$; also, when $a+b>h$, the fused species is $a+b-h$ and the pole of $A_{ab}$ in the strip is at $\theta_1-\theta_2=i\pi(2h-a-b)/h$.
    - *How I will check it:* Carry out the fusion limit symbolically for both $a+b<h$ and $a+b>h$.

20. **Rank 20 — `part3-classical/11-reality.qmd:62` (`#sec-singular-solutions`), `10-solitons.qmd:138` (`#sec-breathers`), `10-solitons.qmd:181` (`#sec-classical-literature`), and `11-reality.qmd:64–75` (`#sec-missing-charges`): Representation of Takács–Watts (`[@tw1999]`) on soliton–breather scattering and excited solitons, and consistency with Part V (`#sec-which-saddles-contribute`) and Part IX (`#sec-critique`, `#sec-sector-status`).**
    - *Claim quoted:* `"In $a_2^{(1)}$, soliton–breather scattering produces complex classical time delays and can drive regular initial data into singular final states [@tw1999]."` and `"The quantum theory also contains excited solitons, which appear as further bound-state poles of the soliton S-matrices. Classical counterparts of some of these are not known [@tw1999] (gap H11 of @sec-critique)."`
    - *Why it is at risk:* Verify that `hep-th/9810006` (`tw1999`) makes both claims as stated, and that the cross-references to `@sec-critique`, `@sec-which-saddles-contribute`, and `@cnj-missing-charge` match what those sections state.
    - *How I will check it:* Read `hep-th/9810006` (`tw1999`) and view the referenced sections in Parts V and IX.

---

## 2. Checked entries (in rank order)

### Rank 1 — `part3-classical/11-reality.qmd:32–41` (`#prp-hirota-energy`), `11-reality.qmd:45–47`, and `10-solitons.qmd:180`

- **Verdict:** **wrong** (the statement of `@prp-hirota-energy` omits two hypotheses that are necessary at complex rapidities and explicitly stated in the cited sources; and `10-solitons.qmd:180` / `11-reality.qmd:45–47` mischaracterize Khastgir and Sasaki's negative-energy solutions).
- **Independent verification:**
  1. *Source check:*
     - Olive, Turok and Underwood (`hep-th/9305160`, Section 8.3, lines 1336–1337 and eq. (8.9)) prove the single-soliton mass formula from the surface term $[\partial_x\ln\tau_j]_{-\infty}^{+\infty}$ and state explicitly: *"In fact this argument holds good even if $e^{\eta_i}$ is complex, providing that its real part is positive"* (i.e. $-\pi/2 < \operatorname{Im}\theta < \pi/2$, or $\operatorname{Re}\cosh\theta > 0$).
     - Takács and Watts (`hep-th/9810006`, Section 4.1, eqs. (4.3)–(4.6), lines 686–698) evaluate the surface term of $\tau = 1 + \sum_i W_i + \dots + (\prod_{i<j} A_{ij})\prod_i W_i$ as $x\to\pm\infty$ under the assumptions that as $x\to+\infty$ every $\exp(m_i x\cosh\theta_i)$ grows ($\operatorname{Re}\cosh\theta_i > 0$) and the top coefficient $\prod_{i<j} A_{ij}(\theta_i-\theta_j)$ is non-zero, so that the single term $(\prod_{i<j}A_{ij})\prod_i W_i$ dominates at $x\to+\infty$ while $1$ dominates at $x\to-\infty$.
     - Khastgir and Sasaki (`hep-th/9507001`, Section 3.2 eq. (3.12), Section 5.2 eqs. (5.12)–(5.13), and Section 6.2 eq. (6.12)): their "complex-energy solitons" (c-solitons) are single solitons with complex rapidity $\theta$ satisfying $\operatorname{Re}\cosh\theta > 0$ (whose real part of energy is positive, $\operatorname{Re}E = M_a\cosh(\operatorname{Re}\theta)\cos(\operatorname{Im}\theta) > 0$, and which blow up at finite time). Their **negative-energy solutions** are *not* Hirota multi-soliton solutions: they are Cauchy-problem solutions with initial data $\vec\phi(x,0) = 0$ and purely imaginary initial velocity $\partial_t\vec\phi(x,0) = i\vec f(x)$ ($\vec f\in L^2(\mathbb R,\mathbb R^r)$), for which $E = -\frac12\int_{-\infty}^\infty \vec f(x)^2\,dx < 0$. Crucially, this initial data satisfies $\Theta\vec\phi(x,t) = -\overline{\vec\phi(x,t)} = \vec\phi(x,t)$ identically! Thus $\Theta$-invariance forces $E$ to be *real*, not positive; positivity of $E$ on $\Theta$-invariant Hirota multi-soliton solutions holds because each constituent has $\operatorname{Re}\cosh\theta_i > 0$.
  2. *Independent computation (`referee_scratch/classical_09_11/check_ch11_reality_breathers.py`, CHECK 2):*
     - For a single soliton of species $a_i$ with $\operatorname{Re}\cosh\theta_i \ne 0$ and no zero of $\tau_j$ on $\mathbb R$, as $x\to+\infty$ we have $E_i\to\infty$ if $\operatorname{Re}\cosh\theta_i > 0$ and $E_i\to 0$ if $\operatorname{Re}\cosh\theta_i < 0$ (while the opposite limit holds as $x\to-\infty$). Consequently,
       $$E = \operatorname{sgn}(\operatorname{Re}\cosh\theta_i)\,M_{a_i}\cosh\theta_i,\qquad P = \operatorname{sgn}(\operatorname{Re}\cosh\theta_i)\,M_{a_i}\sinh\theta_i,$$
       so $\operatorname{Re}E = M_{a_i}\lvert\operatorname{Re}\cosh\theta_i\rvert > 0$ in both regimes! Numerically integrating $T_{00}$ and $-T_{01}$ on $a_2^{(1)}$ ($m=1,\beta=1, M_1=6\sqrt{3}$):
       ```text
       1-soliton Re(cosh th) > 0 (th = 0.3 + 0.5j):
         Numerical int T00 dx       = 9.533598+1.517223j
         Formula M1 cosh(th)        = 9.533598+1.517223j
       1-soliton Re(cosh th) < 0 (th = 0.3 + 2.2j):
         Numerical int T00 dx       = 6.393168-2.558622j
         Formula M1 cosh(th)        = -6.393168+2.558622j
         sgn(Re cosh th) M1 cosh th = 6.393168-2.558622j
       ```
     - When $\prod_{k<l}A_{a_k a_l}(\theta_k-\theta_l) = 0$ — which occurs inside the breather strip of $a_2^{(1)}$ at $u = \pi/6$ ($\theta_1-\theta_2 = i\pi/3$, where $A_{12}(i\pi/3) = \frac{\cos(\pi/3)-\cos(\pi/3)}{\cos(\pi/3)-\cos\pi} = 0$) — the $E_1 E_2$ term in $\tau_j$ vanishes ($\tau_j = 1 + \omega^j E_1 + \omega^{2j}E_2$). Both $E_1$ and $E_2$ grow at the same rate $m_1\cosh\theta_0\cos(\pi/6)$ as $x\to+\infty$, so $[\partial_x\ln(\tau_0\tau_1\tau_2)]_{-\infty}^{+\infty} = 3m_1\cosh\theta_0\cos(\pi/6)$, giving $E = M_1\cosh\theta_0\cos(\pi/6)$, which is **half** of $\sum_{i=1}^2 M_{a_i}\cosh\theta_i = 2M_1\cosh\theta_0\cos(\pi/6)$.
- **Paste-ready replacement text:**
  - In `part3-classical/11-reality.qmd`, lines 32–41, replace:
    ```markdown
    ::: {#prp-hirota-energy}
    **[Established]** For a Hirota solution, the energy and momentum densities are total derivatives [@olive1993a; @olive1993b]. They integrate to

    $$
    E=\sum_iM_{a_i}\cosh\theta_i,\qquad P=\sum_iM_{a_i}\sinh\theta_i,
    $$

    where $(a_i,\theta_i)$ are the species and (possibly complex) rapidities of the constituent solitons [@tw1999, eq. 4.6].
    :::
    ```
    with:
    ```markdown
    ::: {#prp-hirota-energy}
    **[Established]** For a Hirota multi-soliton solution whose constituent solitons $(a_i,\theta_i)$ satisfy $\operatorname{Re}\cosh\theta_i>0$ and $\prod_{k<l}A_{a_ka_l}(\theta_k-\theta_l)\notin\{0,\infty\}$, and for which no $\tau_j$ vanishes on the real line, the energy and momentum densities are total derivatives [@olive1993a; @olive1993b] that integrate to

    $$
    E=\sum_iM_{a_i}\cosh\theta_i,\qquad P=\sum_iM_{a_i}\sinh\theta_i
    $$

    [@olive1993b, sec. 8.3; @tw1999, eq. 4.6]. (If $\operatorname{Re}\cosh\theta_i<0$, the limits $x\to\pm\infty$ exchange roles for $E_i$ and its contribution acquires a minus sign, so $\operatorname{Re}E$ remains positive.)
    :::
    ```
  - In `part3-classical/11-reality.qmd`, lines 47–48, replace:
    ```markdown
    The complex-energy solutions of Khastgir and Sasaki have asymptotic data that are not $\Theta$-invariant.
    ```
    with:
    ```markdown
    The complex-energy single-soliton solutions of Khastgir and Sasaki [@khastgir1996] have complex rapidities whose asymptotic data are not $\Theta$-invariant (while their negative-energy solutions are not Hirota soliton solutions, but $\Theta$-invariant Cauchy-problem configurations with $\vec\phi(x,0)=0$ and purely imaginary $\partial_t\vec\phi(x,0)$).
    ```
  - In `part3-classical/10-solitons.qmd`, line 180, replace:
    ```markdown
    Khastgir and Sasaki [@khastgir1996] found complex-energy and negative-energy solutions and showed that single solitons with complex energy blow up in finite time; @sec-reality-energy-momentum shows that the physical solutions, those invariant under the antilinear symmetry $\Theta$, have real energy.
    ```
    with:
    ```markdown
    Khastgir and Sasaki [@khastgir1996] noted that single solitons at complex rapidity have complex energy and blow up in finite time, and that Cauchy data with $\vec\phi=0$ and purely imaginary $\partial_t\vec\phi$ have negative energy; @sec-reality-energy-momentum shows that Hirota solutions whose constituent data are invariant under the antilinear symmetry $\Theta$ (with $\operatorname{Re}\cosh\theta_i>0$) have real, positive energy.
    ```

---

### Rank 2 — `part3-classical/11-reality.qmd:43–49` (`#sec-reality-energy-momentum`) and `appendices/E-original-results.qmd:8`

- **Verdict:** **wrong** (the $a_2^{(1)}$ spin-2 computation is **correct**, but the general selection rule is stated at `11-reality.qmd:43` without its required additivity/$\xi_i$-independence hypothesis, lacks a statement div/status label/proof in Chapter 11 despite being listed in Appendix E as `**[Theorem, given additivity of the higher charges]**`, and needs the $\Theta$-eigenbasis qualification when an exponent is degenerate).
- **Independent verification:**
  1. *Derivation of the selection rule and its hypotheses:*
     Let $\vec\phi(x,t)$ be a Hirota multi-soliton solution with constituent parameters $\{(a_i,\theta_i,\xi_i)\}_{i=1}^N$. Using MacKay–McGhee (`hep-th/9208057`, Section 3.1–3.4), in every simply-laced affine Toda theory ($a_n^{(1)}$, $d_n^{(1)}$ for both even and odd $n$, and $e_{6,7,8}^{(1)}$), complex conjugation of the $\tau$-functions at each node, $\Theta:\tau_j\mapsto\bar\tau_j$, sends a single soliton with parameters $(a,\theta,\xi)$ to a single soliton of the conjugate species $\bar a$ with parameters $(\bar a,\bar\theta,\bar\xi)$ (in $d_{2m+1}^{(1)}$ the two spinor species have $\delta_{n-1}^{(1)}=-\delta_n^{(1)}=\pm i$, which swap under conjugation; in $e_6^{(1)}$ the two conjugate pairs have coefficients built from $\omega=e^{2\pi i/3}$ and $\omega^2=\bar\omega$, which swap under conjugation). Since $A_{\bar a\bar b}(\bar\theta)=\overline{A_{ab}(\theta)}$, the transformed field $\Theta\vec\phi(x,t) = -\overline{\vec\phi(x,t)}$ is the Hirota multi-soliton solution with parameters $\{(\bar a_i,\bar\theta_i,\bar\xi_i)\}_{i=1}^N$.
     If a local conserved charge $Q=\int_{-\infty}^\infty J_Q(\vec\phi)\,dx$ has density satisfying $J_Q(-\bar{\vec\phi})=\epsilon_Q\overline{J_Q(\vec\phi)}$ with $\epsilon_Q\in\{\pm 1\}$, then
     $$Q[\Theta\vec\phi] = \epsilon_Q \overline{Q[\vec\phi]}.$$
     However, $\Theta$-invariance of the multiset $\{(a_i,\theta_i)\}$ does **not** imply $\Theta\vec\phi = \vec\phi$, because the parameters $\bar\xi_i$ need not match $\xi_{\pi(i)}$ (for an $a_n^{(1)}$ breather, $\xi_2 = \bar\xi_1 + i\pi$, so $\bar\xi_1 = \xi_2 - i\pi \ne \xi_2$). To conclude $Q[\Theta\vec\phi] = Q[\vec\phi]$, one must assume that $Q$ is **additive over the constituents** (or at least depends only on the multiset $\{(a_i,\theta_i)\}$ and is independent of the shifts $\xi_i$), under the regularity hypotheses $\operatorname{Re}\cosh\theta_i > 0$ and $\prod_{k<l}A_{a_k a_l}(\theta_k-\theta_l)\notin\{0,\infty\}$.
  2. *Independent check of the $a_2^{(1)}$ spin-2 current and breather charge (`referee_scratch/classical_09_11/check_ch11_reality_breathers.py`, CHECK 4):*
     In the $\mathbb Z_2$-conjugation basis $\phi_1 = \frac{1}{\sqrt{2}}(-\vec\alpha_0)\cdot\vec\phi$, $\phi_2 = \frac{1}{\sqrt{6}}(\vec\alpha_1-\vec\alpha_2)\cdot\vec\phi$ of $a_2^{(1)}$, we solved $\partial_- J_{+3} + \partial_+ J_{-1} = 0$ symbolically in SymPy:
     $$J_{+3} = \partial_+\phi_2\big((\partial_+\phi_2)^2 - 3(\partial_+\phi_1)^2\big) - \frac{3i}{\sqrt{2}\beta}\big(\partial_+\phi_1\,\partial_+^2\phi_2 - \partial_+^2\phi_1\,\partial_+\phi_2\big),$$
     $$J_{-1} = \frac{3m^2}{2\beta^2}\Big[2e^{i\beta\vec\alpha_0\cdot\vec\phi}\partial_+\phi_2 + e^{i\beta\vec\alpha_1\cdot\vec\phi}(\sqrt{3}\partial_+\phi_1 - \partial_+\phi_2) - e^{i\beta\vec\alpha_2\cdot\vec\phi}(\sqrt{3}\partial_+\phi_1 + \partial_+\phi_2)\Big].$$
     Under $\Theta:\vec\phi\mapsto -\bar{\vec\phi}$, the cubic term in $J_{+3}$ has real coefficients and odd degree $3$, the quadratic derivative term has a purely imaginary coefficient $-3i/(\sqrt{2}\beta)$ and even degree $2$, and $J_{-1}$ has real coefficients times $e^{i\beta\vec\alpha_j\cdot\vec\phi}\partial_+\phi_k$, so every term satisfies $J_{Q_2}(-\bar{\vec\phi}) = -\overline{J_{Q_2}(\vec\phi)}$, confirming $\epsilon_{Q_2} = -1$. Evaluating $Q_2$ numerically on single solitons ($a=1,2$) and on a regular breather ($\theta_0=0.3, u=1.0$ at $t=0$ and $t=0.7$):
     ```text
     th=0.0 [d_- (scales as e^{+2 th})]: Q(a=1) = 38.183766-0.000000j, Q(a=2) = -38.183766+0.000000j
     th=0.4 [d_- (scales as e^{+2 th})]: Q(a=1) = 84.979534-0.000000j, Q(a=2) = -84.979534+0.000000j, Q(a=1)/e^{2th} = 38.183766
     Breather (th0=0.3, u=1.0, t=0.0): Q_2(num) = 0.000000+126.529388j, 2i q_1 e^{2 th0} sin(2u) = 0.000000+126.529388j, diff = 5.00e-14
     Breather (th0=0.3, u=1.0, t=0.7): Q_2(num) = 0.000000+126.529388j, 2i q_1 e^{2 th0} sin(2u) = 0.000000+126.529388j, diff = 6.81e-14
     ```
  3. *Structural check (`11-reality.qmd:43` vs `appendices/E-original-results.qmd:8`):*
     In `appendices/E-original-results.qmd:8`, this result is listed with the status label `**[Theorem, given additivity of the higher charges]**` and cites `(@sec-reality-energy-momentum)`. Per the book's rules (`AGENTS.md`), every item in Appendix E must cite a statement id (`@prp-...` or `@thm-...`) whose status label matches the item's label.
- **Paste-ready replacement text:**
  - In `part3-classical/11-reality.qmd`, lines 43–49, replace:
    ```markdown
    Applying $\Theta$ to this formula gives a selection rule for all conserved charges:

    - If the asymptotic data are $\Theta$-invariant, i.e. the multiset $\{(a_i,\theta_i)\}$ is mapped to itself by $(a,\theta)\mapsto(\bar a,\bar\theta)$, then every conserved charge satisfies $Q=\epsilon_Q\bar Q$, where $J_Q(-\bar{\vec\phi})=\epsilon_Q\overline{J_Q(\vec\phi)}$ defines $\epsilon_Q=\pm1$; in particular $E$ and $P$ are real.
    - A single soliton of a non-self-conjugate species ($a\ne\bar a$) cannot be $\Theta$-invariant alone, because $\Theta$ changes its species. Its energy $M_a\cosh\theta$ is still real when $\theta$ is real, because $M_{\bar a}=M_a$, so the weaker condition $(a,\theta)\mapsto(a,\bar\theta)$ suffices for $E$ and $P$. It does not suffice for charges that distinguish $a$ from $\bar a$. In $a_2^{(1)}$ the spin-2 charge of a single soliton is $q_ae^{2\theta}$, where by Freeman's computation [@freeman1995] the $q_a$ form the eigenvector of the $A_2$ Cartan matrix with eigenvalue $3$, so $q_2=-q_1$, and they can be normalized to be real. If the charge is additive over the constituents, as $E$ and $P$ are, a breather $\{(1,\theta_0+iu),(2,\theta_0-iu)\}$ has spin-2 charge $q_1\big(e^{2(\theta_0+iu)}-e^{2(\theta_0-iu)}\big)=2iq_1e^{2\theta_0}\sin2u$, which is imaginary. This agrees with the selection rule with $\epsilon_Q=-1$. A single soliton of species $1$ at real $\theta$ has real spin-2 charge $q_1e^{2\theta}\ne0$, which is not imaginary: its asymptotic data are not $\Theta$-invariant.
    ```
    with:
    ```markdown
    Applying $\Theta$ gives a selection rule for conserved charges:

    ::: {#prp-theta-selection-rule}
    **[Theorem, given additivity of the higher charges]** Let $Q=\int dx\,J_Q(\vec\phi)$ be a local conserved charge whose density has definite $\Theta$-parity, $J_Q(-\bar{\vec\phi})=\epsilon_Q\overline{J_Q(\vec\phi)}$ with $\epsilon_Q=\pm1$, and assume that on regular Hirota multi-soliton solutions (with $\operatorname{Re}\cosh\theta_i>0$ and $\prod_{k<l}A_{a_ka_l}(\theta_k-\theta_l)\notin\{0,\infty\}$) $Q$ is additive over the constituents, $Q=\sum_i Q_1(a_i,\theta_i)$, and in particular independent of the position parameters $\xi_i$. If the constituent multiset $\{(a_i,\theta_i)\}$ is $\Theta$-invariant, i.e. mapped to itself by $(a,\theta)\mapsto(\bar a,\bar\theta)$, then $Q=\epsilon_Q\bar Q$; in particular $E$ and $P$ ($\epsilon_E=\epsilon_P=+1$) are real and positive.
    :::

    ::: {.proof}
    Under $\Theta:\vec\phi\mapsto-\bar{\vec\phi}$, the $\tau$-functions conjugate at each node, $\tau_j\mapsto\bar\tau_j$, which maps a constituent soliton with parameters $(a_i,\theta_i,\xi_i)$ to one with parameters $(\bar a_i,\bar\theta_i,\bar\xi_i)$. Integrating $J_Q(\Theta\vec\phi)=\epsilon_Q\overline{J_Q(\vec\phi)}$ gives $Q[\Theta\vec\phi]=\epsilon_Q\overline{Q[\vec\phi]}$. Because $Q$ depends only on the multiset $\{(a_i,\theta_i)\}$ and not on the $\xi_i$, $\Theta$-invariance of $\{(a_i,\theta_i)\}$ gives $Q[\Theta\vec\phi]=Q[\vec\phi]$, so $Q=\epsilon_Q\bar Q$. For $E$ and $P$, additivity holds by @prp-hirota-energy and $\operatorname{Re}\cosh\theta_i>0$ gives $E>0$.
    :::

    A single soliton of a non-self-conjugate species ($a\ne\bar a$) cannot be $\Theta$-invariant alone, because $\Theta$ changes its species. Its energy $M_a\cosh\theta$ is still real when $\theta$ is real, because $M_{\bar a}=M_a$, so the weaker condition $(a,\theta)\mapsto(a,\bar\theta)$ suffices for $E$ and $P$. It does not suffice for charges that distinguish $a$ from $\bar a$. In $a_2^{(1)}$ the spin-2 charge of a single soliton is $q_ae^{2\theta}$, where by Freeman's computation [@freeman1995] the $q_a$ form the eigenvector of the $A_2$ Cartan matrix with eigenvalue $3$, so $q_2=-q_1$, and they can be normalized to be real. If the charge is additive over the constituents, as $E$ and $P$ are, a regular breather $\{(1,\theta_0+iu),(2,\theta_0-iu)\}$ ($\frac\pi6<u<\frac\pi2$) has spin-2 charge $q_1\big(e^{2(\theta_0+iu)}-e^{2(\theta_0-iu)}\big)=2iq_1e^{2\theta_0}\sin2u$, which is imaginary. This agrees with @prp-theta-selection-rule with $\epsilon_Q=-1$. A single soliton of species $1$ at real $\theta$ has real spin-2 charge $q_1e^{2\theta}\ne0$, which is not imaginary: its asymptotic data are not $\Theta$-invariant.
    ```
  - In `appendices/E-original-results.qmd`, line 8, replace `(@sec-reality-energy-momentum)` with `(@prp-theta-selection-rule)`.

---

### Rank 3 — `part3-classical/10-solitons.qmd:128–136` (`#sec-breathers`), `11-reality.qmd:28` (`#fig-breather-theta`), and `11-reality.qmd:87` (`#exr-breather-real`)

- **Verdict:** **wrong** (`#exr-breather-real` is **correct**, as is the sine-Gordon case $a=h/2$, but for $a_n^{(1)}$ with $a\ne h/2$, regular breathers built from $(a,\bar a)$ exist only for $\frac{\pi}{2}-\frac{\pi a}{h} < u < \frac{\pi}{2}$ where $A_{a,h-a}(2iu) < 0$, so their classical masses fill $(0, 2M_a\sin\frac{\pi a}{h})$ rather than $(0, 2M_a)$ across $u\in(0,\pi/2)$ as claimed in `10-solitons.qmd:130–136` and `11-reality.qmd:28`; furthermore, `10-solitons.qmd:130–136` omits the type-A breathers of Harder, Iskandar and McGhee [@harder1995] built from $(a,a)$, which are used in `@prp-static-saddles` in Chapter 17).
- **Independent verification:**
  1. *Source check:*
     - Harder, Iskandar and McGhee (`hep-th/9409035`, Section 3.1–3.2, lines 440–612) classify two types of breathers in $a_n^{(1)}$:
       - **Type A** (two solitons of the *same* species $a$ at $\theta_0\pm iu$): $A_{aa}(2iu) = -\frac{\tan^2(\pi a/h)}{\tan^2 u} < 0$ for all $u\in(0,\pi/2)$, giving non-singular oscillating solutions for all $u\in(0,\pi/2)$ with real energy $2M_a\cosh\theta_0\cos u \in (0, 2M_a\cosh\theta_0)$ and topological charge in $\Lambda^{2a\bmod h}\mathbb C^h$ (these are the charged two-soliton bound states used in `@prp-static-saddles` of Chapter 17).
       - **Type B** (conjugate species $a$ and $\bar a = h-a$ at $\theta_0\pm iu$, zero topological charge): with $a\le h/2$, they show (`hep-th/9409035`, lines 525–532) that the breather is non-singular if and only if
         $$\frac{\pi}{2} - \frac{\pi a}{h} < u < \frac{\pi}{2}.$$
     - Takács and Watts (`hep-th/9810006`, Section 4.2, eq. (4.10), lines 714–728) likewise derive for the neutral $(a, h-a)$ breather that non-singularity requires $A_{a,h-a}(2iu) < 0$, i.e. $\cos 2u + \cos(2\pi a/h) < 0$.
  2. *Independent analytic and numerical check (`referee_scratch/classical_09_11/check_ch11_reality_breathers.py`, CHECK 3):*
     From `eq-interaction-coefficient` (`10-solitons.qmd:117`), for $b = h-a$ and $\theta = 2iu$:
     $$A_{a,h-a}(2iu) = \frac{\cosh(2iu) - \cos\frac{\pi(2a-h)}{h}}{\cosh(2iu) - \cos\pi} = \frac{\cos 2u + \cos\frac{2\pi a}{h}}{\cos 2u + 1} = 1 - \frac{\sin^2(\pi a/h)}{\cos^2 u}.$$
     With $\xi_2 = \bar\xi_1 + i\pi$ (`11-reality.qmd:23`), let $\rho = m_a(x\cosh\theta_0 - t\sinh\theta_0)\cos u + \operatorname{Re}\xi_1$ and $\varphi = m_a(x\sinh\theta_0 - t\cosh\theta_0)\sin u + \operatorname{Im}\xi_1$. Then $E_2 = -\bar E_1$ and
     $$\tau_j(x,t) = \big(1 - A_{a,h-a}(2iu)\,e^{2\rho}\big) + 2i\,e^\rho\sin\!\Big(\varphi + \frac{2\pi j a}{h}\Big).$$
     - For $\frac{\pi}{2}-\frac{\pi a}{h} < u < \frac{\pi}{2}$: $A_{a,h-a}(2iu) < 0$, so $\operatorname{Re}\tau_j = 1 + \lvert A_{a,h-a}(2iu)\rvert e^{2\rho} \ge 1 > 0$ everywhere in $(x,t)$. Thus $\tau_j$ never vanishes, has zero winding, and gives a regular neutral breather of rest mass $2M_a\cos u \in \big(0,\,2M_a\sin\frac{\pi a}{h}\big)$.
     - At $u = \frac{\pi}{2}-\frac{\pi a}{h}$: $A_{a,h-a}(2iu) = 0$, the $E_1 E_2$ term disappears, $\tau_j$ vanishes where $e^\rho = \frac12$ and $\sin(\varphi + 2\pi ja/h) = \frac12$, and the integrated energy is $M_a\cos u$, not $2M_a\cos u$.
     - For $0 < u < \frac{\pi}{2}-\frac{\pi a}{h}$: $0 < A_{a,h-a}(2iu) < 1$, so $\operatorname{Re}\tau_j = 0$ on the line $e^{2\rho} = 1/A_{a,h-a}(2iu)$ and $\tau_j$ vanishes whenever $\varphi + 2\pi ja/h \in \pi\mathbb Z$. Even if $\xi_2 - \bar\xi_1$ is varied, any zero-charge breather requires $A_{a,h-a}(2iu)E_1 E_2 > 0$, which for $A_{a,h-a}(2iu) > 0$ forces $E_2 = \bar E_1$ up to a real shift ($\xi_2 = \bar\xi_1$, the case of `@exr-breather-real`), where $\tau_j = 1 + 2e^\rho\cos(\varphi+2\pi ja/h) + A_{a,h-a}(2iu)e^{2\rho}$ is real and attains the negative minimum $1 - 1/A_{a,h-a}(2iu) < 0$ when $\cos(\varphi+2\pi ja/h) = -1$.
     - For `@exr-breather-real` ($\xi_2 = \bar\xi_1$): $\bar\tau_j = \tau_j$ holds for all $u$, and $\tau_j$ changes sign both for $\frac{\pi}{2}-\frac{\pi a}{h} < u < \frac{\pi}{2}$ (where $A<0$, so $\tau_j\to +1$ as $\rho\to-\infty$ and $\tau_j\to -\infty$ as $\rho\to+\infty$) and for $0 < u < \frac{\pi}{2}-\frac{\pi a}{h}$ (where $0<A<1$ and the minimum is $1-1/A<0$). So `@exr-breather-real` is **correct**.
- **Paste-ready replacement text:**
  - In `part3-classical/10-solitons.qmd`, lines 128–136, replace:
    ```markdown
    ## Breathers {#sec-breathers}

    Multi-soliton solutions remain solutions when the rapidities and the parameters $\xi_k$ are continued to complex values. A soliton of species $a$ and one of the conjugate species $\bar a=h-a$, with complex-conjugate rapidities $\theta_0\pm iu$ (so $u$ here is $\frac\pi2$ minus the $u$ of @sec-sg-classical) and suitably related $\xi$'s, form a **breather**: a localized, time-periodic solution of zero topological charge [@hollowood1992; @harder1995]. Its centre of mass moves with rapidity $\theta_0$, its internal oscillation has frequency proportional to $\sin u$, and its mass is

    $$
    M_{\rm breather}=2M_a\cos u,\qquad 0<u<\tfrac\pi2 .
    $$

    In the quantum theory this continuous family is quantized. The lowest breather state of each species $a$ has the mass ratio $\sin\frac{\pi a}h$ of particle $a$, and its S-matrix is the minimal real-coupling S-matrix of Part IV [@hollowood1993; @babujian2013]. In this way the fundamental particles re-emerge inside the soliton theory as quantized breathers, just as in sine-Gordon theory (@sec-sg-breather-particle).
    ```
    with:
    ```markdown
    ## Breathers {#sec-breathers}

    Multi-soliton solutions remain solutions when the rapidities and the parameters $\xi_k$ are continued to complex values. A soliton of species $a\le h/2$ and one of the conjugate species $\bar a=h-a$, with complex-conjugate rapidities $\theta_0\pm iu$ (so $u$ here is $\frac\pi2$ minus the $u$ of @sec-sg-classical) and $\xi_2=\bar\xi_1+i\pi$, form a neutral **breather** (Type B of [@harder1995]): a localized, time-periodic solution of zero topological charge [@hollowood1992; @harder1995; @tw1999]. Here $A_{a,h-a}(2iu)=1-\sin^2(\pi a/h)/\cos^2u$, and regularity ($\tau_j\ne0$ everywhere) requires $A_{a,h-a}(2iu)<0$, which holds for $\frac\pi2-\frac{\pi a}h<u<\frac\pi2$ ( the full interval $(0,\frac\pi2)$ only when $a=h/2$, as in sine-Gordon). Its centre of mass moves with rapidity $\theta_0$, its internal oscillation has frequency proportional to $\sin u$, and its mass is

    $$
    M_{\rm breather}=2M_a\cos u,\qquad \tfrac\pi2-\tfrac{\pi a}h<u<\tfrac\pi2 .
    $$

    (Two solitons of the *same* species $a$ at $\theta_0\pm iu$ have $A_{aa}(2iu)=-\tan^2(\pi a/h)/\tan^2u<0$ and form a charged Type-A breather [@harder1995] of mass $2M_a\cos u$ for all $0<u<\frac\pi2$; @prp-static-saddles uses its static limit.) In the quantum theory the neutral breather family is quantized. The lowest breather state of each species $a$ has the mass ratio $\sin\frac{\pi a}h$ of particle $a$, and its S-matrix is the minimal real-coupling S-matrix of Part IV [@hollowood1993; @babujian2013]. In this way the fundamental particles re-emerge inside the soliton theory as quantized breathers, just as in sine-Gordon theory (@sec-sg-breather-particle).
    ```
  - In `part3-classical/11-reality.qmd`, line 23 and line 28, replace:
    ```markdown
    - **Breathers.** A breather is a two-soliton solution with species $(a,\bar a)$, complex-conjugate rapidities $\theta_1=\theta_0+iu$, $\theta_2=\theta_0-iu$, and $\xi_2=\bar\xi_1+i\pi$ (the $i\pi$ keeps $\tau_j$ away from zero on the real line; see @exr-breather-real).
    ```
    with:
    ```markdown
    - **Breathers.** A neutral breather is a two-soliton solution with species $(a,\bar a)$ ($a\le h/2$), complex-conjugate rapidities $\theta_1=\theta_0+iu$, $\theta_2=\theta_0-iu$ with $\frac\pi2-\frac{\pi a}h<u<\frac\pi2$, and $\xi_2=\bar\xi_1+i\pi$ (where $A_{a,\bar a}(2iu)<0$, so the $i\pi$ keeps $\operatorname{Re}\tau_j\ge1>0$ on the real plane; see @exr-breather-real).
    ```
    and in the caption of `#fig-breather-theta` (`11-reality.qmd:28`), replace:
    ```markdown
    Right: the classical breather mass $2M_a\cos u$, for $u\in(0,\pi/2)$, fills a continuous band.
    ```
    with:
    ```markdown
    Right: for $a=h/2$ (and for Type-A pairs $(a,a)$), the classical breather mass $2M_a\cos u$, for $u\in(0,\pi/2)$, fills $(0,2M_a)$ (for neutral breathers with $a<h/2$, regularity restricts $u$ to $(\frac\pi2-\frac{\pi a}h,\frac\pi2)$ and the mass to $(0,2M_a\sin\frac{\pi a}h)$).
    ```

---

### Rank 4 — `part3-classical/10-solitons.qmd:174` (`#sec-classical-literature`)

- **Verdict:** **correct** (with a minor phrasing note that "and the energies are real" refers to real-rapidity multi-solitons and breathers, as qualified six lines later at line 180).
- **Independent verification:**
  In `referee_scratch/classical_09_11/check_ch11_reality_breathers.py` (CHECK 1), we tested the canonical stress-energy tensor $T_{\mu\nu} = \partial_\mu\vec\phi\cdot\partial_\nu\vec\phi - \eta_{\mu\nu}\mathcal L$ against $(\eta_{\mu\nu}\partial^2 - \partial_\mu\partial_\nu)C$ with $C = -\frac{2}{\beta^2}\sum_j \frac{2}{\vec\alpha_j^2}\ln\tau_j$ symbolically on $a_1^{(1)}$ and numerically on 1-soliton and 2-soliton solutions of $a_2^{(1)}$ at generic complex $(x,t,\theta_i,\xi_i)$:
  ```text
  a_1^(1) 1-soliton: T00 - (-d_x^2 C) = 0, T01 - (-d_t d_x C) = 0, T11 - (-d_t^2 C) = 0
  a_2^(1) 1-soliton (a=1): |T00 - (-C_xx)| = 1.60e-14, |T01 - (-C_xt)| = 1.72e-15, |T11 - (-C_tt)| = 6.90e-15
  a_2^(1) 2-soliton (1,1): |T00 - (-C_xx)| = 1.93e-14, |T01 - (-C_xt)| = 1.57e-15, |T11 - (-C_tt)| = 1.41e-14
  a_2^(1) 2-soliton (1,2): |T00 - (-C_xx)| = 1.21e-14, |T01 - (-C_xt)| = 4.89e-15, |T11 - (-C_tt)| = 2.46e-14
  ```
  This agrees with Olive, Turok and Underwood (`olive1993a; hep-th/9305160`) and Aratyn et al. (`hep-th/9212086`, eqs. (2.21)–(2.25)): on every vertex-operator soliton solution, the conformal stress tensor vanishes identically, so the canonical stress tensor equals the improvement divergence $(\eta_{\mu\nu}\partial^2 - \partial_\mu\partial_\nu)C$ pointwise.

---

### Rank 5 — `part3-classical/10-solitons.qmd:22–25` (`#sec-imaginary-coupling-vacua`), `10-solitons.qmd:104` (`#sec-soliton-masses`), and `10-solitons.qmd:158–160` (`#exr-lie-duality`)

- **Verdict:** **correct**.
- **Independent verification:**
  1. *Source check:* Kneipp and Olive (`hep-th/9404030`, Section 6, lines 951–986) and Olive, Turok and Underwood (`hep-th/9305160`, Section 7, lines 1190–1204) show that for twisted affine algebras folded by a centre/Weyl automorphism, single parent solitons survive with $M_{\langle a\rangle} = M_a \propto m_a$, whereas for untwisted non-simply-laced affine algebras folded by a Dynkin diagram automorphism fixing node $0$, the folded soliton of orbit $O_a$ is a coincident orbit of $\lvert O_a\rvert = 2/\vec\alpha_a^2$ parent solitons, giving $M_a = \frac{2}{\vec\alpha_a^2}\frac{2h}{\beta^2}m_a$, whose ratios are the particle mass ratios of $(g^\vee)^{(1)}$.
  2. *Independent computation (`referee_scratch/classical_09_11/check_ch10_solitons.py`, CHECK 6):*
     We computed the eigenvalues of $\mathsf M^2 = m^2\sum_j n_j \vec\alpha_j\vec\alpha_j^{\mathsf T}$ and the folded soliton masses for $b_3^{(1)}$, $c_3^{(1)}$, $g_2^{(1)}$, and $f_4^{(1)}$:
     ```text
     b_3^(1) particle masses: [1.414214 1.414214 2.44949 ]
     c_3^(1) soliton mass ratios (scaled by sqrt(2)): [1.414214 1.414214 2.44949 ]
     c_3^(1) particle masses: [1.       1.732051 2.      ]
     b_3^(1) soliton mass ratios (scaled by 1/sqrt(2)): [1.       1.732051 2.      ]
     g_2^(1) particle masses: [1.414214 2.44949 ] ratio m_2/m_1 = 1.7320508075688772
     g_2^(1) soliton mass ratio M_2/M_1 = 1.7320508075688776 (matches sqrt(3)!)
     f_4^(1) particle masses: [1.126033 1.59245  2.175328 3.076378]
     f_4^(1) soliton mass ratios: [1.59245  2.252065 3.076378 4.350655]
     ```
     For $a_{2n}^{(2)}$ folded from $a_{2n}^{(1)}$, every soliton is a coincident pair ($\lvert O_a\rvert = 2$) of species $(a, 2n+1-a)$, so every soliton mass is multiplied by the same factor $2$ and remains proportional to $m_a = 2m\sin\frac{\pi a}{2n+1}$.

---

### Rank 6 — `part3-classical/10-solitons.qmd:76–77` (`#sec-vertex-solitons`)

- **Verdict:** **wrong** (the $c_2^{(1)}$ static $\tau$-functions are **correct**, as verified symbolically, and $n_j^{\rm parent} = n_j^\vee$ holds for the eight non-simply-laced foldings of `@tbl-foldings`, but (a) for the self-dual twisted algebra $a_{2n}^{(2)}$ folded from $a_{2n}^{(1)}$ (`06-affine-algebras.qmd:113`, `18-jordan-chains.qmd:350`), the degree of a coincident pair ($\lvert O\rvert = 2$) is $\lvert O\rvert n_j^{\rm parent} = 2$ at every node whereas $\lvert O\rvert n_j^\vee$ gives $4$ at nodes $0,\dots,n-1$ and $2$ at node $n$; and (b) writing $\tau_1 = \tau_3$ in $c_2^{(1)}$ uses the parent $a_3^{(1)}$ node numbering $0,1,2,3$ without explanation).
- **Independent verification:**
  1. *Symbolic check of $c_2^{(1)}$ (`referee_scratch/classical_09_11/check_ch10_solitons.py`, CHECK 2):*
     With $E = \exp(\sqrt{2}mx)$, $\tau_0 = 1 + 2E + \frac12 E^2$, $\tau_1 = 1 + \frac12 E^2$, $\tau_2 = 1 - 2E + \frac12 E^2$, we substituted into `eq-hirota-bilinear` for $c_2^{(1)}$ ($n_0=1, n_1=2, n_2=1$, $\vec\alpha_0^2=\vec\alpha_2^2=2, \vec\alpha_1^2=1$, so $n_j\vec\alpha_j^2/2 = 1$ for all $j$, with $-a_{1,0}=-a_{1,2}=2$ and $-a_{0,1}=-a_{2,1}=1$):
     ```text
     c_2^(1) static soliton Hirota residuals: j=0 -> 0, j=1 -> 0, j=2 -> 0
     ```
  2. *Comparison of $n_j^{\rm parent}$ and $n_j^\vee$ across foldings:*
     In `@hyp-D` (`18-jordan-chains.qmd:323`), the degree of $\tau_j$ for a folded orbit $O$ of parent solitons in a simply-laced parent $P^{(1)}$ is $\lvert O\rvert n_j^{\rm parent}$, where $n_j^{\rm parent}$ is the parent Kac label at node $j$. For the 4 untwisted foldings ($b_n^{(1)}, c_n^{(1)}, f_4^{(1)}, g_2^{(1)}$) and 4 twisted foldings ($a_{2n-1}^{(2)}, d_{n+1}^{(2)}, e_6^{(2)}, d_4^{(3)}$) in `@tbl-foldings`, each node orbit $O_j$ consists of mutually orthogonal parent roots ($\vec\alpha_j^2 = 2/\lvert O_j\rvert$) and $n_j = \lvert O_j\rvert n_j^{\rm parent}/k$, so $n_j^\vee = k n_j \vec\alpha_j^2 / 2 = n_j^{\rm parent}$ holds identically (and matches the level of $L(\Lambda_j)$ in Olive–Turok–Underwood `hep-th/9305160`, Section 7, lines 1200–1203). However, for $a_{2n}^{(2)}$ folded from $a_{2n}^{(1)}$ (`06-affine-algebras.qmd:113`), the orbit $\{n, n+1\}$ has adjacent roots so $\vec\alpha_n^2 = 1/2$, $n_j = (1,2,\dots,2)$, and $n_j^\vee = 2n_j\vec\alpha_j^2/2 = (2,2,\dots,2,1) \ne n_j^{\rm parent} = (1,1,\dots,1)$.
- **Paste-ready replacement text:**
  In `part3-classical/10-solitons.qmd`, lines 76–77, replace:
  ```markdown
  For folded theories the relevant label is $n_j^\vee$, the level of $L(\Lambda_j)$. A coincident orbit $O$ of $\lvert O\rvert$ parent solitons has degree $\lvert O\rvert n_j^\vee$. In $c_n^{(1)}$ every $n_j^\vee$ is 1, so a coincident pair of species $a$ and $2n-a$ has degree 2 at every node, including the long end nodes, whose Kac label is 1, and the species-$n$ soliton has degree 1 at every node, including the short nodes, whose Kac label is 2. For $c_2^{(1)}$, folded from $a_3^{(1)}$, the static pair of species 1 and 3 has $\tau_0=1+2E+\frac12E^2$, $\tau_1=\tau_3=1+\frac12E^2$ and $\tau_2=1-2E+\frac12E^2$.
  ```
  with:
  ```markdown
  For theories folded from a simply-laced parent (@tbl-foldings), the relevant label at node $j$ is the parent Kac label $n_j^{\rm parent}$, which for every untwisted and twisted pair in @tbl-foldings equals $n_j^\vee$, the level of $L(\Lambda_j)$ (for $a_{2n}^{(2)}$ folded from $a_{2n}^{(1)}$, every $n_j^{\rm parent}=1$). A coincident orbit $O$ of $\lvert O\rvert$ parent solitons has degree $\lvert O\rvert n_j^{\rm parent}$. In $c_n^{(1)}$ every $n_j^\vee=n_j^{\rm parent}$ is 1, so a coincident pair of species $a$ and $2n-a$ has degree 2 at every node, including the long end nodes, whose Kac label is 1, and the species-$n$ soliton has degree 1 at every node, including the short nodes, whose Kac label is 2. For $c_2^{(1)}$, folded from $a_3^{(1)}$, the static pair of species 1 and 3 has $\tau_0=1+2E+\frac12E^2$, $\tau_1=1+\frac12E^2$ (from $\tau_1=\tau_3$ in $a_3^{(1)}$) and $\tau_2=1-2E+\frac12E^2$.
  ```

---

### Rank 7 — `part3-classical/11-reality.qmd:58–62` (`#sec-singular-solutions`) and `11-reality.qmd:90–92` (`#exr-singular-points`)

- **Verdict:** **wrong** (the double pole with zero residue is **correct**, but (a) a single soliton with real $\theta$ depends only on $x\cosh\theta - t\sinh\theta$, so its real zeros are never isolated points in $\mathbb R^2$; (b) a double pole on the real line is not Lebesgue/Cauchy-PV integrable, so local conserved charges are unaffected only when the integration contour is deformed into $\mathbb C$ around the pole; and (c) the **topological charge** is *not* unaffected when a zero of $\tau_j$ crosses the real $x$-axis, as it jumps by the coroot $\pm\vec\alpha_j^\vee$).
- **Independent verification (`referee_scratch/classical_09_11/check_ch11_reality_breathers.py`, CHECK 5):**
  1. *Isolated points vs curves:* If $\tau_j(x,t) = 0$ at $(x_0,t_0)\in\mathbb R^2$, writing $\tau_j = \operatorname{Re}\tau_j + i\operatorname{Im}\tau_j = 0$ gives two real equations in two real variables $(x,t)$, whose solutions are generically isolated points by the inverse function theorem whenever the Jacobian $\operatorname{Im}(\overline{\partial_x\tau_j}\,\partial_t\tau_j) \ne 0$ (which holds e.g. for a single soliton at complex rapidity $\operatorname{Im}\theta\notin\pi\mathbb Z$, or a generic multi-soliton). For a single soliton at *real* rapidity $\theta$, $\tau_j$ is a function of the single real variable $x\cosh\theta - t\sinh\theta$, so $\operatorname{Im}(\overline{\partial_x\tau_j}\,\partial_t\tau_j) = 0$ identically and any real zero forms a straight line in $\mathbb R^2$.
  2. *Double pole with zero residue:* Near a zero $x_0$ of $\tau_j(\cdot,t)$ of multiplicity $m_0 \ge 1$, $\tau_j(x,t) = (x-x_0)^{m_0} g(x)$ with $g(x_0)\ne 0$. Then
     $$\partial_x\ln\tau_j = \frac{m_0}{x-x_0} + \frac{g'(x)}{g(x)},\qquad \partial_x^2\ln\tau_j = -\frac{m_0}{(x-x_0)^2} + \partial_x\left(\frac{g'(x)}{g(x)}\right).$$
     Because $\partial_x^2\ln\tau_j$ is the $x$-derivative of $\partial_x\ln\tau_j$, the coefficient of $(x-x_0)^{-1}$ in $T_{00} = \frac{2}{\beta^2}\sum_j \frac{2}{\vec\alpha_j^2}\partial_x^2\ln\tau_j$ is **identically zero**. Thus an integration contour in the complex $x$-plane bypassing $x_0$ above or below gives the same integral $[\partial_x\ln\tau_j]_{-\infty}^{+\infty}$ (whereas the real-line integral $\int_{\mathbb R} T_{00}\,dx$ diverges at $x_0$).
  3. *Jump of the topological charge:* Since $\vec Q_{\rm top} = \frac{\beta}{2\pi}\big(\vec\phi(+\infty)-\vec\phi(-\infty)\big) = -\frac{1}{2\pi}\sum_j \vec\alpha_j^\vee [\operatorname{Im}\ln\tau_j]_{-\infty}^{+\infty}$, when a simple zero $x_0(t)$ of $\tau_j$ crosses the real $x$-axis as $t$ passes $t_0$, the winding angle $[\operatorname{Im}\ln\tau_j]_{-\infty}^{+\infty}$ jumps by $\pm 2\pi$, so $\vec Q_{\rm top}$ jumps by $\mp\vec\alpha_j^\vee$!
- **Paste-ready replacement text:**
  - In `part3-classical/11-reality.qmd`, lines 60–61, replace:
    ```markdown
    Zeros of $\tau_j$ in the real $(x,t)$ plane are generically isolated points: one complex condition, two real conditions. Non-generic data, such as a single soliton with tuned $\operatorname{Im}\xi$ or a breather with $\xi_2=\bar\xi_1$, give curves of zeros. There the energy density has double poles with zero residue, so integrated charges are unaffected, and the Hirota data continue analytically through the singularity.
    ```
    with:
    ```markdown
    Zeros of $\tau_j$ in the real $(x,t)$ plane are generically isolated points (one complex condition, two real conditions, when $x$ and $t$ enter with linearly independent complex coefficients, as at complex rapidity or in multi-soliton collisions). Non-generic data, such as a real-rapidity single soliton with tuned $\operatorname{Im}\xi$ or a breather with $\xi_2=\bar\xi_1$, give curves of zeros. At a zero of $\tau_j$, the energy density has a double pole in $x$ with zero residue, so local conserved charges defined by deforming the integration contour around the pole are unaffected, and the Hirota data continue analytically through the singularity (whereas the topological charge defined along the real axis jumps by a coroot $\pm\vec\alpha_j^\vee$ when a zero crosses the real line).
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

---

### Rank 8 — `part3-classical/09-lagrangian.qmd:90–96` (`#eq-spin3-current`) and `09-lagrangian.qmd:168–170` (`#exr-spin3`)

- **Verdict:** **wrong** (`eq-spin3-current` is **correct** for sinh-Gordon, but (a) at line 96, `"For sine-Gordon, replace $\beta^2$ by $-\beta^2$"` is false unless $\cosh(\beta\phi)$ is also replaced by $\cos(\beta\phi)$, i.e. replace $\beta$ by $i\beta$; and (b) with $\partial_\pm = \partial_t \pm \partial_x$, $\partial_+$ scales on a right-moving state of rapidity $\theta$ as $e^{-\theta}$, so the charge $Q_3 = \int dx\,(T_4 + \mathcal X_2)$ built from $\partial_+$ in `eq-spin3-current` has spin $-3$ and scales as $e^{-3\theta}$, evaluating on the sine-Gordon kink to $-\frac{16m^3}{3\beta^2}e^{-3\theta}$, whereas the charge scaling as $e^{+3\theta}$ is built from $\partial_-$).
- **Independent verification (`referee_scratch/classical_09_11/check_ch09_lagrangian.py`, CHECK 1, and `check_rank10_exponents.py`):**
  1. *Conservation check:*
     ```text
     1. Sinh-Gordon conservation residual d_- T_4 + d_+ X_2 = 0
     2. Sine-Gordon with ONLY beta^2 -> -beta^2 (keeping cosh): d_- T_4 + d_+ X_2 = m**2*p_p*(beta*p_p**2*sin(beta*phi) + beta*p_p**2*sinh(beta*phi) - 2*p_pp*cos(beta*phi) + 2*p_pp*cosh(beta*phi)) != 0
     3. Sine-Gordon with beta -> i*beta (beta^2 -> -beta^2 AND cosh -> cos): d_- T_4 + d_+ X_2 = 0
     ```
  2. *Lorentz spin and evaluation on the sine-Gordon kink:*
     For the sine-Gordon kink $\phi(x,t) = \frac{4}{\beta}\arctan\exp(u)$ with $u = m(x\cosh\theta - t\sinh\theta)$ (`eq-sg-kink` in `03-sine-gordon.qmd:16`), $\partial_+ u = (\partial_t + \partial_x)u = m(\cosh\theta - \sinh\theta) = m e^{-\theta}$. With $\partial_u\phi = \frac{2}{\beta}\operatorname{sech}u$, $\partial_u^2\phi = -\frac{2}{\beta}\operatorname{sech}u\tanh u$, and $\cos(\beta\phi) = 1 - 2\operatorname{sech}^2 u$:
     $$T_4 = (\partial_+^2\phi)^2 - \frac{\beta^2}{4}(\partial_+\phi)^4 = \frac{4m^4}{\beta^2}e^{-4\theta}\big(\operatorname{sech}^2 u - 2\operatorname{sech}^4 u\big),$$
     $$\mathcal X_2 = m^2\cos(\beta\phi)(\partial_+\phi)^2 = \frac{4m^4}{\beta^2}e^{-2\theta}\big(\operatorname{sech}^2 u - 2\operatorname{sech}^4 u\big).$$
     Since $\int_{-\infty}^\infty (\operatorname{sech}^2 u - 2\operatorname{sech}^4 u)\,du = 2 - \frac83 = -\frac23$ and $dx = du/(m\cosh\theta)$:
     $$Q_{-3} = \int_{-\infty}^\infty (T_4 + \mathcal X_2)\,dx = -\frac{8m^3}{3\beta^2\cosh\theta}\big(e^{-4\theta} + e^{-2\theta}\big) = -\frac{16m^3}{3\beta^2}e^{-3\theta}.$$
     (Swapping $\partial_+ \leftrightarrow \partial_-$ gives the spin-$+3$ current with $Q_{+3} = -\frac{16m^3}{3\beta^2}e^{+3\theta}$; if $\mathcal X_2$ is shifted by $-m^2(\partial_\pm\phi)^2$ so that $\mathcal X_2 = m^2(\cos\beta\phi - 1)(\partial_\pm\phi)^2$, then $Q_{\pm 3} = -\frac{32m^3}{3\beta^2}e^{\pm 3\theta}$.)
- **Paste-ready replacement text:**
  In `part3-classical/09-lagrangian.qmd`, line 96, replace:
  ```markdown
  The charge $Q_3=\int dx\,\big(T_4+\mathcal X_2\big)$ built from this current has spin 3, and on a one-particle state it takes the value $\chi\,m^3e^{3\theta}$, as in (-@eq-charge-eigenvalue). For sine-Gordon, replace $\beta^2$ by $-\beta^2$.
  ```
  with:
  ```markdown
  Because $\partial_+=\partial_t+\partial_x$ multiplies a mode of rapidity $\theta$ by $e^{-\theta}$, the charge $\int dx\,\big(T_4+\mathcal X_2\big)$ built from (-@eq-spin3-current) has spin $-3$ and scales as $e^{-3\theta}$; exchanging $\partial_+\leftrightarrow\partial_-$ gives the spin-$3$ charge $Q_3$, which on a one-particle state takes the value $\chi\,m^3e^{3\theta}$, as in (-@eq-charge-eigenvalue). For sine-Gordon, replace $\beta$ by $i\beta$, so $\beta^2\to-\beta^2$ and $\cosh(\beta\phi)\to\cos(\beta\phi)$.
  ```

---

### Rank 9 — `part3-classical/10-solitons.qmd:82–90` (`#sec-topological-charges`), `10-solitons.qmd:92` (`#fig-topological-charges`), and `10-solitons.qmd:154–156` (`#exr-charges-a3`)

- **Verdict:** **correct**.
- **Independent verification (`referee_scratch/classical_09_11/check_ch10_solitons.py`, CHECK 1):**
  Writing the simple roots of $a_n^{(1)}$ ($h=n+1$) in the orthonormal basis $\vec e_1,\dots,\vec e_h$ of $\mathbb R^h$ as $\vec\alpha_j = \vec e_j - \vec e_{j+1}$ ($1\le j\le n$) and $\vec\alpha_0 = \vec e_h - \vec e_1$, for any $c_j = -\frac{1}{2\pi}[\vartheta_j]_{(-\pi,\pi)}$ with $\vartheta_j = \frac{2\pi j a}{h} + \operatorname{Im}\xi$ we have
  $$\frac{\beta}{2\pi}\vec\phi(+\infty) = \sum_{j=0}^{h-1} c_j\vec\alpha_j = \sum_{k=1}^h (c_k - c_{k-1})\vec e_k = \sum_{k=1}^h \Big(I_k - \frac{a}{h}\Big)\vec e_k,$$
  where $I_k \in \{0,1\}$ counts whether $\vartheta_j$ crosses an odd multiple of $\pi$ between $k-1$ and $k$, and $\sum_{k=1}^h I_k = a$. Thus $\frac{\beta}{2\pi}\vec\phi(+\infty)$ has $a$ components equal to $1 - a/h$ and $h-a$ components equal to $-a/h$, which is a weight of $\Lambda^a\mathbb C^h$ (and not $\Lambda^{h-a}\mathbb C^h$). As $\operatorname{Im}\xi$ increases by $2\pi$, the reduced phases jump at $h/\gcd(a,h)$ values, $\gcd(a,h)$ of them at once, cycling through $h/\gcd(a,h)$ weights of $\Lambda^a\mathbb C^h$. All weights for $a_2^{(1)}, a_3^{(1)}, a_4^{(1)}$ were verified in `check_ch10_solitons.py`.

---

### Rank 10 — `part3-classical/09-lagrangian.qmd:84–88` (`#sec-conserved-charges`)

- **Verdict:** **correct**.
- **Independent verification (`referee_scratch/classical_09_11/check_rank10_exponents.py`):**
  With $h = \sum_j n_j$ (sum of coprime Kac labels), the period $kh$ and exponents modulo $kh$ match Kac (`@kac1990`, Chapter 14 and Tables Aff 1–3) for every untwisted and twisted affine algebra:
  - $a_2^{(2)}$: $k=2, h=3, kh=6$, exponents $1, 5 \pmod 6$;
  - $d_4^{(3)}$: $k=3, h=4, kh=12$, exponents $1, 5, 7, 11 \pmod{12}$, matching $g_2^{(1)}$ ($1, 5 \pmod 6$);
  - $e_6^{(2)}$: $k=2, h=9, kh=18$, exponents $1, 5, 7, 11, 13, 17 \pmod{18}$, matching $f_4^{(1)}$ ($1, 5, 7, 11 \pmod{12}$);
  - $a_{2n}^{(2)}$: $k=2, h=2n+1, kh=4n+2$, odd integers not divisible by $2n+1$;
  - $b_n^{(1)}, c_n^{(1)}, a_{2n-1}^{(2)}, d_{n+1}^{(2)}$: all odd integers (verified in `check_rank10_exponents.py` that $E_1^s \in \mathfrak{so}(2N)$ is non-zero and $\sigma$-invariant for every odd $s \ge 1$, including $s = 2n-1$ in $d_{2n}^{(1)} \to a_{2n-1}^{(2)}$).

---

### Rank 11 — `part3-classical/09-lagrangian.qmd:122–138` (`#sec-masses-couplings`, `#eq-area-rule`, `#exr-shg-cubic`)

- **Verdict:** **wrong** (`#eq-area-rule` and `#exr-shg-cubic` are **correct**, but `09-lagrangian.qmd:130` states that for degenerate masses in $(a_n, d_n, e_6)$ the charge eigenvectors $\vec e_a$ are complex and that *"In a real orthonormal basis of a degenerate pair neither (-@eq-area-rule) nor @prp-dorey-rule holds"*; this is true for complex-conjugate pairs with $\bar a \ne a$, i.e. $a_n^{(1)}$ ($a\ne h/2$), $d_{2m+1}^{(1)}$ spinors, and $e_6^{(1)}$, but **false** for $d_{2m}^{(1)}$ with even $n=2m$, where the centre $\mathbb Z_2\times\mathbb Z_2$ has real characters $\pm 1$, every particle is self-conjugate $\bar a = a$, and the charge eigenbasis of the spinor pair $(s,s')$ *is* a real orthonormal basis in which both `eq-area-rule` and `@prp-dorey-rule` hold).
- **Independent verification (`referee_scratch/classical_09_11/check_ch09_lagrangian.py`, CHECK 4):**
  In $d_4^{(1)}$ ($h=6$, even $n=4$), the $\mathbb Z_2\times\mathbb Z_2$ centre characters are real ($\pm 1$) and the character eigenbasis is the real orthonormal basis $\vec u_1 = (1,0,0,0)$, $\vec u_s = (0,0,0,1)$, $\vec u_{s'} = (0,1,1,0)/\sqrt{2}$, $\vec u_2 = (0,1,-1,0)/\sqrt{2}$. Computing $C_{abc}$ and $\frac{4\beta}{\sqrt{6}}\Delta_{abc}$ ($m=1,\beta=1$) in this real orthonormal basis gives:
  ```text
  d_4^(1) (h=6) in the REAL orthonormal Z_2 x Z_2 character eigenbasis:
    |C_{1,1,s'}| = 1.414214,  (4/sqrt(6))*Delta = 1.414214,  diff = 8.88e-16
    |C_{1,1,2(central)}| = 1.414214,  (4/sqrt(6))*Delta = 1.414214,  diff = 2.00e-15
    |C_{s,s,s'}| = 1.414214,  (4/sqrt(6))*Delta = 1.414214,  diff = 8.88e-16
    |C_{s,s,2(central)}| = 1.414214,  (4/sqrt(6))*Delta = 1.414214,  diff = 2.00e-15
    |C_{s',s',2(central)}| = 1.414214,  (4/sqrt(6))*Delta = 1.414214,  diff = 2.22e-15
    |C_{2(central),2(central),2(central)}| = 4.242641,  (4/sqrt(6))*Delta = 4.242641,  diff = 8.88e-16
  ```
- **Paste-ready replacement text:**
  In `part3-classical/09-lagrangian.qmd`, line 130, replace:
  ```markdown
  For degenerate masses ($a_n$, $d_n$, $e_6$) the $\vec e_a$ must be charge eigenvectors, the complex eigenvectors of the centre of the affine diagram (of $\mathbb Z_h$ for $a_n$), normalized by $\vec e_a\cdot\vec e_{\bar a}=1$ so that $\phi_{\bar a}=\bar\phi_a$. Particle $a$ is the eigenvector with the centre character of node $a$. In a real orthonormal basis of a degenerate pair neither (-@eq-area-rule) nor @prp-dorey-rule holds.
  ```
  with:
  ```markdown
  For degenerate masses ($a_n$, $d_n$, $e_6$) the $\vec e_a$ must be charge eigenvectors, the eigenvectors of the centre of the affine diagram (of $\mathbb Z_h$ for $a_n$; see @sec-hirota-form-dn1-en1), normalized by $\vec e_a\cdot\vec e_{\bar a}=1$ so that $\phi_{\bar a}=\bar\phi_a$ (which are complex for conjugate pairs $a\ne\bar a$ in $a_n$, $d_{2m+1}$ and $e_6$, and real orthonormal for even $d_{2m}$, where $\bar a=a$). Particle $a$ is the eigenvector with the centre character of node $a$. In a real orthonormal basis of a complex-conjugate pair ($a\ne\bar a$), neither (-@eq-area-rule) nor @prp-dorey-rule holds.
  ```

---

### Rank 12 — `part3-classical/09-lagrangian.qmd:142–160` (`#sec-dorey-rule`, `#prp-dorey-rule`, `#fig-dorey-d4`, `#exr-dorey-a2`)

- **Verdict:** **correct**.
- **Independent verification (`referee_scratch/classical_09_11/check_ch09_lagrangian.py`, CHECK 4):**
  For white node $a$, $s_{\rm black}\vec\lambda_a = \vec\lambda_a$ and $s_{\rm white}\vec\lambda_a = \vec\lambda_a - \vec\alpha_a$, so $\vec\gamma_a = (1 - s_{\rm white}s_{\rm black})\vec\lambda_a = \vec\alpha_a$; for black node $a$, $s_{\rm white}\vec\lambda_a = \vec\lambda_a$ and $s_{\rm black}\vec\lambda_a = \vec\lambda_a - \vec\alpha_a$, so $\vec\gamma_a = \vec\lambda_a - s_{\rm white}(\vec\lambda_a - \vec\alpha_a) = s_{\rm white}\vec\alpha_a$. In $D_4$ ($h=6$), the four Coxeter orbits $\Omega_1, \Omega_{2(c)}, \Omega_{3(s)}, \Omega_{4(s')}$ have 6 roots each, and searching all triples of orbits confirms that $\Omega_a + \Omega_b + \Omega_c \ni 0$ holds for $(1,1,2)$, $(s,s,2)$, $(s',s',2)$, $(1,s,s')$, and $(2,2,2)$, matching the non-zero $C_{abc}$ in CHECK 4 and the triangle $\vec\gamma_1 + \vec\gamma_1' + \vec\gamma_c = 0$ with side ratios $1:1:\sqrt{3}$ in `@fig-dorey-d4`.

---

### Rank 13 — `part3-classical/10-solitons.qmd:122` (`#sec-multi-solitons`)

- **Verdict:** **correct** (with a minor notation observation in Section 4 below that $E_1$ denotes both the exponential of soliton 1 and the energy of soliton 1 in the same paragraph).
- **Independent verification (`referee_scratch/classical_09_11/check_ch10_solitons.py`, CHECK 3):**
  For $\theta_1 > \theta_2$, along soliton 1's trajectory $x(t) = (\tanh\theta_1)t + x_0$, the exponent of $E_2$ is $m_{a_2}\cosh\theta_2(\tanh\theta_1 - \tanh\theta_2)t + \text{const}$, which goes to $-\infty$ as $t\to-\infty$ ($E_2\to 0$) and to $+\infty$ as $t\to+\infty$ ($E_2\to\infty$). Thus soliton 1's position shifts by $\Delta x = -\frac{\ln A_{a_1a_2}(\theta_1-\theta_2)}{m_{a_1}\cosh\theta_1} > 0$, giving a time advance $\Delta t = -\Delta x/v_1 < 0$. With $E_1 = M_{a_1}\cosh\theta_1$, $p_1 = M_{a_1}\sinh\theta_1$, and $M_{a_1} = \frac{2h}{\beta^2}m_{a_1}$, we get $E_1\Delta x = -p_1\Delta t = -\frac{2h}{\beta^2}\ln A_{a_1a_2}(\theta_1-\theta_2)$. In the rest frame of soliton 2 ($\theta_2=0, \theta=\theta_1$), $dE_1/d\theta = p_1$, so $\frac{d\delta_{ab}}{d\theta} = p_1\Delta t = \frac{2h}{\beta^2}\ln A_{ab}(\theta)$. For sine-Gordon ($h=2$, $\beta_{\rm SG}^2 = 2\beta^2$), $A_{11}(\theta) = \frac{\cosh\theta-1}{\cosh\theta+1} = \tanh^2\frac{\theta}{2}$, so $\frac{2h}{\beta^2}\ln A_{11}(\theta) = \frac{4}{\beta^2}\ln\tanh^2\frac{\theta}{2} = \frac{16}{\beta_{\rm SG}^2}\ln\tanh\frac{\theta}{2}$, matching the small-$\xi$ limit ($\pi/\xi \simeq 8\pi/\beta_{\rm SG}^2$) of $-i\frac{d}{d\theta}\ln S_0(\theta)$ from `eq-sg-s0` (`03-sine-gordon.qmd:70`).

---

### Rank 14 — `part3-classical/09-lagrangian.qmd:114–117` (`#sec-masses-couplings`) and `#exr-dn-masses` (`09-lagrangian.qmd:172–174`)

- **Verdict:** **correct**.
- **Independent verification (`referee_scratch/classical_09_11/check_ch09_lagrangian.py`, CHECK 3):**
  We computed the eigenvalues of $\mathsf M^2 = m^2\sum_j n_j \vec\alpha_j\vec\alpha_j^{\mathsf T}$ for $d_4^{(1)}, d_5^{(1)}, d_6^{(1)}$ and $e_8^{(1)}$:
  ```text
  d_4^(1) (h=6): masses and sum rule 2h=12 match to 4.44e-16
  d_5^(1) (h=8): masses and sum rule 2h=16 match to 2.22e-16
  d_6^(1) (h=10): masses and sum rule 2h=20 match to 4.44e-16
  e_8^(1) (h=30): max diff between eigenvalues and table formulas = 8.88e-16, sum m_a^2 = 60.000000 (2h = 60)
  ```

---

### Rank 15 — `part3-classical/10-solitons.qmd:10`

- **Verdict:** **correct**.
- **Independent verification (`referee_scratch/classical_09_11/check_ch10_solitons.py`, CHECK 4):**
  In $a_{2n-1}^{(1)}$ ($h=2n$), for species $a=n$ we have $\omega^{jn} = e^{2\pi i j n/(2n)} = (-1)^j$. At $\operatorname{Im}\xi = \pi/2$, $E = iR$ with $R = e^{\operatorname{Re}\xi + m_n(x\cosh\theta - t\sinh\theta)} > 0$, so $\tau_{2k} = 1 + iR$ and $\tau_{2k+1} = 1 - iR = \bar\tau_{2k}$. Using $\sum_{k=0}^{n-1}\vec\alpha_{2k+1} = -\sum_{k=0}^{n-1}\vec\alpha_{2k}$, we obtain
  $$\vec\phi(x,t) = \frac{i}{\beta}\sum_{j=0}^{2n-1}\vec\alpha_j\ln\tau_j = -\frac{2}{\beta}\left(\sum_{k=0}^{n-1}\vec\alpha_{2k}\right)\arctan R,$$
  which is strictly real and proportional to the fixed vector $\sum_{k=0}^{n-1}\vec\alpha_{2k}$ times a sine-Gordon kink.

---

### Rank 16 — `part3-classical/09-lagrangian.qmd:57–79` (`#sec-lax-pair`, `#eq-lax-pair`, `#eq-zero-curvature`, `#exr-lax`)

- **Verdict:** **correct**.
- **Independent verification (`referee_scratch/classical_09_11/check_ch09_lagrangian.py`, CHECK 2):**
  Using $[E_{\vec\alpha_i}, E_{-\vec\alpha_j}] = \delta_{ij}\frac{2\vec\alpha_j}{\vec\alpha_j^2}\cdot\vec H$ and $[\vec H, E_{\pm\vec\alpha_j}] = \pm\vec\alpha_j E_{\pm\vec\alpha_j}$, the $E_{\pm\vec\alpha_j}$ terms in $[\partial_+ + \mathcal A_+, \partial_- + \mathcal A_-]$ cancel identically, and the $\vec H$ coefficient is $-\big(\beta\partial_+\partial_-\vec\phi + m^2\sum_j n_j\vec\alpha_j e^{\beta\vec\alpha_j\cdot\vec\phi}\big)\cdot\vec H = 0$. Under $\vec\phi_\vartheta(x^+,x^-) = \vec\phi(e^{-\vartheta}x^+, e^{\vartheta}x^-)$, $\partial_\pm\vec\phi_\vartheta(x) = e^{\mp\vartheta}(\partial_\pm\vec\phi)(x')$, so $\mathcal A_\pm[\vec\phi_\vartheta](x;\lambda) = e^{\mp\vartheta}\mathcal A_\pm[\vec\phi](x'; e^\vartheta\lambda)$. The $2\times 2$ matrix check for $a_1^{(1)}$ (`@exr-lax`) gives `F_curv[0,0] = -sqrt(2)*beta*phi_pm/2 - 2*m^2*sinh(sqrt(2)*beta*phi) = 0`, recovering the sinh-Gordon equation.

---

### Rank 17 — `part3-classical/10-solitons.qmd:68–79` (`#sec-vertex-solitons`) vs `part2-algebra/06-affine-algebras.qmd:165–171`

- **Verdict:** **correct**.
- **Independent verification:**
  In `06-affine-algebras.qmd:168`, the matrix element is written with $Q_k W_k(x,t) F^{a_k}(z_k)$, while `10-solitons.qmd:72` writes $Q_k F^{a_k}(z_k)$ with the spacetime dependence absorbed into $F^{a_k}(z_k) = e^{-E_1 x^+ + E_{-1}x^-}\hat F^{a_k}(z_k)e^{E_1 x^+ - E_{-1}x^-}$ (or into $Q_k$), as in Olive, Turok and Underwood (`hep-th/9305160`, Section 6, eq. (6.3)). In `hep-th/9305160`, Section 7, lines 1135–1168, the central charge of $\hat F^a(z)$ in the principal Heisenberg algebra is $1$ on a level-$1$ representation, so on $L(\Lambda_j)$ of level $x = n_j$ (for simply-laced $\hat g$) the normal-ordering zero prevents more than $x = n_j$ factors of $\hat F^a(z)$ at the same rapidity from being non-zero, giving $(F^a(z))^{n_j+1} = 0$ on $L(\Lambda_j)$.

---

### Rank 18 — `part3-classical/10-solitons.qmd:27` (`#sec-imaginary-coupling-vacua`)

- **Verdict:** **correct**.
- **Independent verification:**
  Checked against `part2-algebra/07-quantum-groups.qmd` and `part6-soliton-smatrices/20-quantum-group-symmetry.qmd` (`eq-q-physical`). In Chapter 20 (`eq-q-physical`), Part VI uses the opposite coproduct to Chapter 7 and fixes $q_{\rm VI} = e^{-4\pi^2 i/\beta^2} = -e^{-i\pi\omega}$, which is $q^{-1}$ for $q = e^{4\pi^2 i/\beta^2}$ of `10-solitons.qmd:27` (for twisted algebras $\hat g = X_n^{(k)}$, $\beta^2$ is replaced by $\beta_{\rm lit}^2 = \beta^2/k$ as stated in `06-affine-algebras.qmd:153`).

---

### Rank 19 — `part3-classical/10-solitons.qmd:126` (`#sec-multi-solitons`) and `10-solitons.qmd:162–164` (`#exr-fusion`)

- **Verdict:** **wrong** (the trigonometric identity $m_a e^{i\pi b/h} + m_b e^{-i\pi a/h} = m_{a+b}$ is **correct** for $a+b < h$, but (a) both `10-solitons.qmd:126` and `@exr-fusion` omit the condition $a+b < h$: when $a+b > h$, such as $a=b=2$ in $a_2^{(1)}$ where $2+2\equiv 1\pmod 3$, $2m\sin\frac{\pi(a+b)}{h} < 0$ and the fused soliton of species $a+b-h$ arises instead at $\theta_1 = \theta - i\pi(h-b)/h$, $\theta_2 = \theta + i\pi(h-a)/h$ where $\theta_1-\theta_2 = i\pi(2h-a-b)/h$; and (b) the limit must be taken by sending $\theta_1-\theta_2 \to i\pi(a+b)/h$ jointly with $\xi_1,\xi_2\to-\infty$ so that $A_{ab}(\theta_1-\theta_2)e^{\xi_1+\xi_2} \to e^\xi$, rather than fixing $\theta_1-\theta_2 = i\pi(a+b)/h$ first where $A_{ab}=\infty$).
- **Independent verification (`referee_scratch/classical_09_11/check_ch10_solitons.py`, CHECK 5):**
  For $1\le a,b\le h-1$:
  $$2m\sin\frac{\pi a}{h}e^{i\pi b/h} + 2m\sin\frac{\pi b}{h}e^{-i\pi a/h} = 2m\sin\frac{\pi(a+b)}{h},$$
  which equals $+m_{a+b} > 0$ when $a+b < h$, vanishes when $a+b = h$, and equals $-m_{a+b-h} < 0$ when $a+b > h$. For $a+b > h$, the fused species is $c = a+b-h \in \{1,\dots,h-2\}$, and replacing $(a,b)$ by the conjugate shifts $-(h-b)$ and $+(h-a)$ gives
  $$2m\sin\frac{\pi a}{h}e^{-i\pi(h-b)/h} + 2m\sin\frac{\pi b}{h}e^{i\pi(h-a)/h} = 2m\sin\frac{\pi(a+b-h)}{h} = +m_{a+b-h} > 0$$
  at $\theta_1 - \theta_2 = -i\pi(2h-a-b)/h$ (where $\cosh(\theta_1-\theta_2) = \cos\frac{\pi(a+b)}{h}$ so the denominator of $A_{ab}$ vanishes).
- **Paste-ready replacement text:**
  - In `part3-classical/10-solitons.qmd`, line 126, replace:
    ```markdown
    **Fusion.** The two-soliton solution also contains the single solitons as residues. At $\theta_1-\theta_2=i\pi(a+b)/h$, the denominator of $A_{ab}$ vanishes. Take $\theta_1=\theta+i\pi b/h$ and $\theta_2=\theta-i\pi a/h$. Then $E_1E_2$ has the exponent of a single soliton of species $a+b$ with rapidity $\theta$, because $m_ae^{i\pi b/h}+m_be^{-i\pi a/h}=m_{a+b}$, and the phases combine to $\omega^{j(a+b)}$. Letting $\xi_1,\xi_2\to-\infty$ while keeping $A_{ab}e^{\xi_1+\xi_2}$ fixed reduces (-@eq-n-soliton) to the single soliton of species $a+b$ (@exr-fusion).
    ```
    with:
    ```markdown
    **Fusion.** The two-soliton solution also contains the single solitons as residues. For $a+b<h$, the denominator of $A_{ab}$ vanishes at $\theta_1-\theta_2=i\pi(a+b)/h$. Taking $\theta_1\to\theta+i\pi b/h$ and $\theta_2\to\theta-i\pi a/h$, $E_1E_2$ acquires the exponent of a single soliton of species $a+b$ with rapidity $\theta$, because $m_ae^{i\pi b/h}+m_be^{-i\pi a/h}=m_{a+b}$, and the phases combine to $\omega^{j(a+b)}$. Sending $\xi_1,\xi_2\to-\infty$ jointly with $\theta_1-\theta_2\to i\pi(a+b)/h$ so that $A_{ab}(\theta_1-\theta_2)e^{\xi_1+\xi_2}\to e^\xi$ reduces (-@eq-n-soliton) to the single soliton of species $a+b$ (@exr-fusion; for $a+b>h$, take $\theta_1\to\theta-i\pi(h-b)/h$ and $\theta_2\to\theta+i\pi(h-a)/h$ to obtain species $a+b-h$).
    ```
  - In `part3-classical/10-solitons.qmd`, lines 162–164 (`#exr-fusion`), replace:
    ```markdown
    ::: {#exr-fusion}
    Verify the fusion identity $m_ae^{i\pi b/h}+m_be^{-i\pi a/h}=m_{a+b}$ for $a_n^{(1)}$ and carry out the limit that reduces the two-soliton $\tau_j$ (-@eq-n-soliton) to the single soliton of species $a+b$.
    :::
    ```
    with:
    ```markdown
    ::: {#exr-fusion}
    For $a_n^{(1)}$ with $a+b<h$, verify the fusion identity $m_ae^{i\pi b/h}+m_be^{-i\pi a/h}=m_{a+b}$ and carry out the joint limit $\theta_1-\theta_2\to i\pi(a+b)/h$, $\xi_1,\xi_2\to-\infty$ with $A_{ab}(\theta_1-\theta_2)e^{\xi_1+\xi_2}\to e^\xi$ that reduces the two-soliton $\tau_j$ (-@eq-n-soliton) to the single soliton of species $a+b$.
    :::
    ```

---

### Rank 20 — `part3-classical/11-reality.qmd:62` (`#sec-singular-solutions`), `10-solitons.qmd:138` (`#sec-breathers`), `10-solitons.qmd:181` (`#sec-classical-literature`), and `11-reality.qmd:64–75` (`#sec-missing-charges`)

- **Verdict:** **correct**.
- **Independent verification:**
  1. *Source check (`hep-th/9810006`, Takács and Watts 1999):*
     - In Section 4.3–4.4 (lines 800–960), Takács and Watts compute classical soliton–breather scattering in $a_2^{(1)}$ and show that the classical time delay of the soliton is complex, shifting $\operatorname{Im}\xi$ so that a regular incoming soliton + breather state is generically driven into a singular outgoing state where $\tau_j$ has real zeros.
     - In Section 3.5 and Section 5 (lines 590–660, 1080–1160), they compare the excited-soliton bound-state poles of the Hollowood quantum soliton S-matrix with the semiclassical quantization of the classical breather / two-soliton solutions and conclude that several quantum excited solitons lack known classical counterparts.
  2. *Cross-reference check:*
     Checked `part5-semiclassics/17-thimbles.qmd` (`#sec-which-saddles-contribute`, `@cnj-missing-charge`, `@prp-static-saddles`) and `part9-foundations/44-twelve-gaps.qmd` (`#sec-critique`, `#sec-sector-status`): all cross-references, status labels (`[Open]`, `[Conjecture]`), and descriptions match Part III accurately.

---

## 3. Unchecked entries and why

None. All 20 entries in Section 1 were checked in full.

---

## 4. Other observations (not verified)

1. **`part3-classical/10-solitons.qmd:122` (`#sec-multi-solitons`): Overload of the symbol $E_1$ in a single paragraph.**
   In the first sentence of the "Time delays" paragraph, $E_1$ and $E_2$ denote the Hirota exponentials of solitons 1 and 2 (`"either $E_2\to0$, and $\tau_j\simeq1+\omega^{ja_1}E_1$, or $E_2\to\infty$"`). Five sentences later in the same paragraph, $E_1$ denotes the energy of soliton 1 (`"the shift is $E_1\Delta x=-p_1\Delta t=-(2h/\beta^2)\ln A_{a_1a_2}$, with $E_1,p_1$ the energy and momentum of soliton 1"`), despite the note at line 61 (`"In this chapter and the next... $E$ is the exponential in (-@eq-an-soliton), not an energy"`). Writing $\mathcal E_1$ or $M_{a_1}\cosh\theta_1\,\Delta x$ for the energy of soliton 1 at line 122 would avoid the collision.
2. **`part3-classical/09-lagrangian.qmd:90` (`#eq-spin3-current`) vs `09-lagrangian.qmd:165` (`#exr-lax`): Different normalizations of $\beta$ and $m$ for sinh-Gordon within the same chapter.**
   At line 90, sinh-Gordon is written with $V = \frac{m^2}{\beta^2}(\cosh\beta\phi - 1)$ and $\partial_+\partial_-\phi = -\frac{m^2}{\beta}\sinh\beta\phi$ (the standalone convention of Chapter 3, where the particle mass is $m$), whereas in `#eq-imaginary-potential` (line 24) and `@exr-lax` (line 165), $a_1^{(1)}$ has $\vec\alpha_1 = -\vec\alpha_0 = \sqrt{2}$, $V_{\rm real} = \frac{2m^2}{\beta^2}(\cosh(\sqrt{2}\beta\phi)-1)$, and particle mass $m_1 = 2m$. Adding a parenthetical `(in the convention of @sec-sg-classical, with $\beta_{\rm SG}$ written as $\beta$)` at line 90 would prevent confusion when comparing `@exr-lax` and `@exr-spin3`.
