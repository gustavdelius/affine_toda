# Draft report: chapter 39, "Solitons and particles through a defect"

Drafted 2026-10-12 by Claude (Opus) from `briefs/draft_ch39.md` plus the later instructions (reuse ch. 38's `eq-dc-ar-delay`, derive the type-I₊ delays, ch. 38 obligations). Files: `part8-defects/39-defect-solitons.qmd` (236 lines, 38.5 kB; ch. 31 is 36.7 kB), `textbook_code/defect_solitons_an.py` (9 checks, all PASS, about 10 s). No other file was edited and no git command that changes the index or the working tree was run.

## 1. Chapter summary

Sections (all skeleton ids kept):

- **Intro + conventions.** Notation $\gamma_a=\pi a/h$, $\hat M=2hm/\beta^2$, $M_a=2\hat M\sin\gamma_a$; complex rapidities in the strip $|\operatorname{Im}\theta|<\pi/2$ (a species-$a$ soliton at $\theta+i\pi$ is the reversed $\bar a$ soliton).
- **`sec-an-soliton-defect`.** The *class* $p\in\Lambda_W/\Lambda_R\cong\mathbb Z_h$ of a jump; new `eq-ds-class-shift`: a jump of class $p$ multiplies all $f_j$ by $\omega^{-p}$, so it acts as $\theta_d\to\theta_d\pm2\pi ip/h$ (types I$_\mp$). `eq-an-defect-delay` is the general statement: $z^{\rm I_-}_{a,p}=z_a(\theta-\theta_d-2\pi ip/h)$, $z^{\rm I_+}_{a,p}=z_a(\theta_d-2\pi ip/h-\theta)$, with $z_a$ = ch. 38's `eq-dc-ar-delay` (cited, plus its sinh form). `prp-an-defect-delay` **[Established; numerically checked for $a_2$ and $a_3$, all species, both types and all classes]**: the four items (dictionary, $|z|>1$, limits $1$ and $\omega^a$ = rotation $\vec\psi=\pi\vec\phi$, Bäcklund). Proof: order-$E$ = linearized problem at $\theta-i\pi/2$; classes by the shift; I₊ by time reversal. Remark: $\theta_d\to-\infty$ turns I₋ into the invertible line of the diagram symmetry $\pi$ (a_r version of the NOT gate and of `prp-defect-twist-limit` (3)).
- **`sec-defect-topological-charge`.** New `eq-ds-class-energy` (stored $(\mathcal B_d,P_d)$ of class $p$; with CZ constants on the mass shell $\hat M$; real iff $2p\equiv0$ mod $h$). New `eq-ds-delay-phase` ($\arg z_a=\arg(e^{-2(\theta-\theta_d)}+e^{2i\gamma_a})$, book form of cz2007 (2.22)). `prp-stored-charge` **[Established; numerically checked]**: bookkeeping (jump changes by $\lambda_{\rm in}-\lambda_{\rm out}\in\Lambda_R$), cycle = $\pi$-orbit with boundaries $2\pi g/h$, selection rule $0\le k\le a/g$ steps ($\bar a/g$ backwards for $a>h/2$), both bounds attained, no missing charges created. Remark: the charge changes only through the singular (curve-of-zeros) soliton.
- **`sec-defect-absorption`.** New `eq-ds-absorption` $\theta_{\rm abs}=\theta_d\mp i(\pi/2-\pi(a+2p)/h)$. New `prp-ds-absorption` **[Theorem]**: absorption/emission configurations; envelope velocity $\tanh\theta_d$ (direction rule for all species and classes); orientation by the sign of $\sin(\pi(a+2p)/h)$; stored $M_a(\cosh,\sinh)\theta_{\rm abs}$; real iff $a+2p\equiv h/2$ mod $2h$ (species $h/2-2q$ in class $q$; e.g. species 1 of $a_5$ in class 1). The a_r analogue of "stores a soliton of rapidity $\theta_d$": class $a$ stores species $a$ at $\theta_d-i(\pi/2-\gamma_a)$.
- **`sec-defect-multi`.** Multi-solitons ($E_k\to z_kE_k$, $A_{ab}$ unchanged), classical factorization (triangle relation), two defects multiply, new `eq-ds-bootstrap` = cz2007 (2.25) in book form.
- **`sec-defect-reality`.** Θ maps class $p\to-p$ and conjugates delays; breathers stay Θ-invariant with real displacement; charged classes $p\not\equiv0,h/2$ are complex sectors (advanced solitons if $\cos(2\pi p/h)<0$, non-unimodular linear factors); odd sector for $h$ even = class $h/2$ = $\sigma\to-\sigma$, $z_a(u+i\pi)=1/z_{\bar a}(u)$; for $h$ odd $\sigma\to-\sigma$ is a different defect.
- **`sec-classical-transmission-factors`.** `eq-an-linear-transmission`: $R=0$, $\mathbb T^{\rm cl}_k=(e^{-u}+e^{i\gamma_k})/(e^{-u}+e^{-i\gamma_k})=z_k(u+i\pi/2)$ (I₋), $z_k(\theta_d-\theta-i\pi/2)$ (I₊). `prp-an-linear-transmission` **[Numerical]** (5 items incl. crossing $\bar{\mathbb T}_k(\theta)=\mathbb T_{\bar k}(i\pi-\theta)$, limits, class-$p$ factors, cz2007 (7.5)/(7.11)). Pre-registered frozen-soliton test, then `prp-defect-frozen-soliton` **[Theorem]**: (1) same-type pairs reproduce $\kappa A_{cb}^{\pm1}$ only for $b=h/2$; (2) a type-I₋ defect at $\vartheta=\theta_c+i(\pi/2+\gamma_c)$ and a type-I₊ at $\bar\vartheta$ give $\omega^bA_{cb}(\theta-\theta_c)$ and $\omega^bX_{cb}$ for every $b$, and carry $M_c(\cosh,\sinh)\theta_c$ with CZ constants; (3) real-σ realization by classes $4p\equiv h+2c$ mod $4h$ when $h+2c\equiv0$ mod 4. Outcome reported as "a qualified yes, in a twisted form" (soliton = I₋ ⊗ I₊ Bianchi factors up to the $\pi$-line).
- **Exercises (8):** `exr-ds-class`, `exr-ds-linear`, `exr-ds-time-reversal`, `exr-ds-selection`, `exr-ds-absorption`, `exr-ds-bootstrap`, `exr-ds-breather`, `exr-ds-frozen`.
- **Notes** with history and four corrections (below).

**New ids** (all prefixed `ds-`, only in this chapter): `eq-ds-class-shift`, `eq-ds-class-energy`, `eq-ds-delay-phase`, `eq-ds-absorption`, `eq-ds-bootstrap`, `prp-ds-absorption`, `exr-ds-class`, `exr-ds-linear`, `exr-ds-time-reversal`, `exr-ds-selection`, `exr-ds-absorption`, `exr-ds-bootstrap`, `exr-ds-breather`, `exr-ds-frozen`. All cross-references resolve to existing ids (grep), all citation keys exist in `references.bib`; no new bib keys needed.

## 2. Script output (`python3 textbook_code/defect_solitons_an.py`, about 10 s)

```
PASS delays (a_2, a_3, all species, both types, every class p): ... at 7 times [res 1e-13; class-0 delay in class p != 0: >= 3.98]; least-squares z agrees to 5e-16; time reversal I_-(sigma) -> I_+(1/sigma) [2e-14]; B_d(phi, psi + 2 pi lambda_p/beta; sigma) = B_d(phi, psi; sigma e^(-+2 pi i p/h)) + const, root shifts invariant [3e-14]; = BCZ (6.2) at theta -> -theta = cz2007 (2.20) at a -> h - a [2e-15]
PASS stored charge: static coweight jumps solve the conditions; class p stores ...; real iff 2p = 0 mod h (a_2..a_5, both types); ... jump change (t-continuation at x = 0) = lambda_in - lambda_out [2e-15], jumps in the weight lattice; arg z_a monotone ...; sector boundaries 2 pi gcd(a,h)/h apart; max steps along the cycle = a/gcd(a,h) [((3,1),(1,1)), ((3,2),(1,1)), ((4,1),(1,1)), ((4,2),(1,1)), ((4,3),(1,1)), ((5,1),(1,1)), ((5,2),(2,2)), ((5,3),(2,2)), ((5,4),(1,1))]; theta >> theta_d: lambda_out = pi lambda_in
PASS cz2007 (2.22): the sector boundaries are 2 gamma_a apart only when a divides h (e.g. a_4, a = 2: spacing 2 pi/5, not 4 pi/5; a species-2 soliton can move two steps)
PASS absorption/emission at theta_abs ...: (soliton, const) and (const, soliton) solve the conditions [res 4e-14]; ... [rel 6e-16]; ... E + B_d, P + P_d constant by quadrature through the real absorptions a_3 (a = 2, p = 0), a_5 (a = 1, p = 1) [rel 6e-10]; theta - theta_d -> +inf: z_a -> omega^a, psi = pi phi (the rotation line)
PASS several solitons: ... [res 1e-13]; two defects (I_- then I_+) multiply the delays [1e-14]; z_a(u + i pi b/h) z_b(u - i pi a/h) = z_(a+b)(u) (cz2007 (2.25)), h = 3..8 [3e-15]
PASS reality: ... breather [res 4e-13]; Theta-image of a class-p solution solves the conditions (class -p) [5e-14]; ...; real rapidities: delayed (|z| > 1) iff cos(2 pi p/h) > 0, advanced otherwise
PASS linear waves (real coupling; a_1..a_4, both types, all particles): R = 0, no mixing, |T| = 1 [3e-15]; closed forms [3e-15]; crossing/unitarity [5e-15]; a_1 [3e-16]; limits [4e-13]; T(theta_d), poles, zeros [9e-16]; Hessian vs finite differences [5e-11]
PASS cz2007 in their conventions (D^T, rho_k, m = 1): corrected (7.11) [1e-15], = book T_(h-k), = (7.5) for a_2 [1e-15]; the printed (7.11) misses by >= 0.50; about a class-p defect (imaginary coupling) T_k at theta_d + 2 pi i p/h [3e-10], unimodular only if 2p = 0 mod h
PASS frozen soliton: X_ab(i sinh u) = A_ab(u + i pi/2) and T_b(u) = z_b(u + i pi/2) [9e-16]; same type iff b = h/2 (zero/pole matching, h = 2..8; a_2 grid search: best deviation 1.64); mixed pair [9e-14], energy with BCZ constants [1e-15]; real-sigma classes (a_1 c = 1, a_3 c = 2, a_5 c = 1, 3) [3e-16], on fields [8e-14]
   runtime 10 s
```

No check failed. One development note: a first version of check (c) assumed the class always goes $p\to p+a$; it failed for $\sin(\pi(a+2p)/h)<0$, which led to the orientation statement in `prp-ds-absorption` item 2 (the object is then the reversed $\bar a$ soliton and the classes go $p+a\to p$). The grid search originally returned `inf` because the $u$-grid contained the zero of $A_{11}$; fixed.

## 3. Verification record

**By script** (`defect_solitons_an.py`, check numbers as printed):
1. `eq-an-defect-delay`, `prp-an-defect-delay`, `eq-ds-class-shift`: delays for all species of $a_2$, $a_3$, both types, all classes, real and complex $\theta$; control (class-0 delay fails in class $p\ne0$); least-squares fits (I₊ class 1, I₊ class 0 of $a_3$, I₋ class 2 of $a_3$); time reversal; coweight/root shifts of $\mathcal B_d$ ($r=2,3,4$); dictionary to bcz2004b (6.2), cz2007 (2.20).
2. `eq-ds-class-energy`, `prp-stored-charge`, `eq-ds-delay-phase`: static class energies vs cz2009a (4.8), mass shell, reality; charges by $x$-continuation vs formula; jump change by $t$-continuation; weight-lattice jumps; arg monotonicity/range; boundary spacing; max steps $=a/g$ attained ($a_2,a_3,a_4$); $\theta\gg\theta_d$ gives $\pi\lambda$; cz2007 spacing claim fails for $(h,a)=(5,2)$.
3. `prp-ds-absorption`, `eq-ds-absorption`: configurations at $\theta_{\rm abs}$ ($a_2,a_3,a_5$; both types; $\theta_d=\pm0.6$); orientation sign; $\Delta\mathcal B_d$, $\Delta P_d$; class values before/after; zero of class-$(p+a)$ delay; envelope velocity; real-rapidity rule for $h\le8$; quadrature of $E+\mathcal B_d$, $P+P_d$ (Simpson) through real absorptions; rotation limit $\vec\psi=\pi\vec\phi$.
4. `sec-defect-multi`, `eq-ds-bootstrap`: two-soliton passage; two defects in a row; bootstrap identity $h=3..8$.
5. `sec-defect-reality`: Θ on delays (all classes); breather; Θ-image in class $-p$; conjugate class energies; $|z|>1$ iff $\cos(2\pi p/h)>0$.
6.–7. `eq-an-linear-transmission`, `prp-an-linear-transmission`: full mode-space linear solve (both sides) for $r=1..4$, both types; analytic Hessian vs finite differences; values, poles, zeros; limits; cz2007 (7.11) corrected vs printed and (7.5), in their conventions; class-$p$ linear factors (finite-difference Hessian at imaginary coupling).
8. `prp-defect-frozen-soliton`: continuation identities; zero/pole matching for $h\le8$; grid search ($a_2$); mixed-pair identities incl. $X$; energy; real-σ class realization ($a_1,a_3,a_5$) and the two-defect field configuration ($a_3$, $a_5$).

**By hand (drafting agent):** $\vec\alpha_j\cdot\mathsf D\vec\lambda_p$ values (also printed numerically in a scratch run); the order-$E$ argument that $z_a(\theta)=\mathbb T^{\rm cl}_a(\theta-i\pi/2)$; $|z_a|^2$ formula; the energy identity for $\theta_{\rm abs}$; the zero/pole proof of `prp-defect-frozen-soliton`; the gradient formulas of `exr-ds-linear` (also checked numerically in a scratch run, $10^{-15}$); $z_a(u+i\pi)=1/z_{\bar a}(u)$; the statement that a fused I₋I₊ pair at the CZ normalization carries $M_c$ but $2\hat M+M_c$ in the book's normalization.

**Unchecked / not derived here:**
- The general $N$-soliton passage (only $N=2$ checked; the text says the general case follows from the Bäcklund property and Bianchi permutability, "which we have not written out for $a_r$").
- That the selection-rule transition happens through the singular (curve-of-zeros) soliton: follows from the charge formula, not scripted.
- The reading of charged classes as "complex sectors excluded by Θ" in the quantum theory (an analogy, phrased as such).
- Whether the fused I₋I₊ pair is one of the type-II defects of `sec-type-two-atft`, and whether its quantum transmission matrix contains the soliton S-matrix (stated as open).
- cz2007 Sec. 6 tunnelling discussion and cz2009a (5.8) (quoted for the quantum matrices; not re-derived).

## 4. Corrections to the sources (in the Notes)

1. cz2007 (7.11): denominator $ie^{-\theta}m_k+\sigma(1-\omega^k)$ (confirms the digest), and the formula assumes $m=1$ (otherwise $\sigma\to m\sigma$). Checked in their own conventions.
2. cz2007 Sec. 2 after (2.22): sector boundaries are $2\pi\gcd(a,h)/h$ apart, not $2\gamma_a$ (agree only if $a\,|\,h$), as cz2009a Sec. 4 already notes; quantified: up to $a/\gcd(a,h)$ steps (two for species 2 of $a_4$).
3. cz2007 after (2.23) and cz2009a after (4.10): "real-rapidity absorption only for self-conjugate solitons" holds for root-lattice jumps only; class-$p$ defects absorb species $a$ at real rapidity iff $a+2p\equiv h/2$ mod $2h$.
4. cz2009a after (4.8): static energies are also real when $h$ is even and the classes differ by $h/2$.

(cz2007 (2.21), $-\coth$ not $+\coth$, is ch. 38's correction and is not repeated; `prp-an-defect-delay` item 1 uses $-\coth$ consistently.)

## 5. Requests

**Edits elsewhere (integration phase):**
- `textbook_code/README.md`: row for `defect_solitons_an.py` (ch. 39; ids `eq-an-defect-delay`, `prp-an-defect-delay`, `eq-ds-class-shift`, `eq-ds-class-energy`, `prp-stored-charge`, `prp-ds-absorption`, `eq-ds-bootstrap`, `eq-an-linear-transmission`, `prp-an-linear-transmission`, `prp-defect-frozen-soliton`; 9 checks; about 10 s).
- `part8-defects/38-defect-conditions.qmd`, end of `sec-defect-backlund-an`: it promises that "@sec-defect-solitons develops ... the absence of absorption at real rapidity except for self-conjugate solitons". With charged jumps this is only true in class 0; suggest "except for self-conjugate solitons, unless the defect carries a charged jump (@prp-ds-absorption)". Also its sentence "the phase of $z_a$ is what moves the topological charge ... to a neighbouring weight" could say "along its cycle of charges (@prp-stored-charge)".
- `part8-defects/37-defect-sine-gordon.qmd` (`prp-half-soliton` discussion): optional pointer to `prp-defect-frozen-soliton` (the a_r version is an I₋I₊ pair up to the π-line).
- `part3-classical/10-solitons.qmd` `sec-topological-charges`: optional pointer to `prp-stored-charge` (the defect moves charges along the π-orbit and never reaches missing charges).
- `part5-semiclassics/16-dhn.qmd` after `eq-an-transmission`: optional pointer to `eq-an-linear-transmission` ($\mathbb T^{\rm cl}=z$ continued, as $X=A$ continued).
- `part8-defects/36-defect-scattering.qmd` `sec-topological-lines` / `prp-twist-defect`: optional sentence that the $\theta_d\to-\infty$ limit of the classical a_r defect is the invertible line of the diagram rotation $\pi$.

**Conventions brief (`part8_defects_conventions.md`) additions:**
- Jump class $p\in\mathbb Z_h$; `eq-ds-class-shift`: class $p$ ≡ $\theta_d\to\theta_d\pm2\pi ip/h$ (I$_\mp$), i.e. $\sigma\to\sigma\omega^{\mp p}$; $h$ even, $p=h/2$ is $\sigma\to-\sigma$ (the a_r odd sector).
- Type-I₊ a_r delays $z^{\rm I_+}_{a,p}(\theta)=z_a(\theta_d-2\pi ip/h-\theta)$ (time reversal), linear factors $\mathbb T^{\rm I_+}_k=z_k(\theta_d-\theta-i\pi/2)$.
- $\mathbb T^{\rm cl}_k(\theta)=z_k(\theta-\theta_d+i\pi/2)$ for I₋ (book mode $k$: $\pi\rho_k=\omega^k\rho_k$); cz2007's $T_k$ is book $\mathbb T_{h-k}$.
- Strip convention: complex soliton rapidities in $|\operatorname{Im}\theta|<\pi/2$; outside it the soliton is the reversed conjugate species.
- Class energies `eq-ds-class-energy`; CZ normalization subtracts $\hat M(\cosh,\sinh)\theta_d$.
- Notation of ch. 39: $\gamma_a=\pi a/h$ (not a coupling), $\hat M=2hm/\beta^2$, frequency written $\nu$ (since $\omega=e^{2\pi i/h}$).

**Obligations for later chapters:**
- **Ch. 41** (`sec-an-soliton-transmission`, `sec-transmission-semiclassics`): (i) the a_r transmission matrices should respect the selection rule of `prp-stored-charge` in the classical limit (at most $a/g$ steps, direction of $\pi$); compare cz2009a (5.8) and the cz2007 (5.9) "tunnelling" asymmetry, with book species labels ($a\to h-a$); (ii) the Wigner/semiclassical test with $\ln|z_a|$ and the phase $\arg z_a$ (`eq-ds-delay-phase`), for both types; (iii) the class-$p$ sectors: the quantum matrices in charged sectors should be the class-0 ones at $\theta_d\pm2\pi ip/h$ (the a_r generalization of ch. 37's $s=1$ intertwiner at $\theta_d-i\pi$); (iv) test whether a fused I₋I₊ pair at $\vartheta,\bar\vartheta$ reproduces the soliton S-matrix times the $\pi$-line (the a_r analogue of cz2010 Sec. 4), and whether it is a type-II defect; (v) the I₋/I₊ assignment to $U_q(\mathfrak b^\mp)$ for a_r with the classical $\theta_d\to-\infty$ limit = $\pi$-line.
- **Ch. 41** (`sec-an-defect-bound-states`): the classical poles `eq-ds-absorption`, including the real absorptions by charged classes ($a+2p\equiv h/2$ mod $2h$), against the cz2009a (6.13)–(6.15) unstable-state windows (whose "soliton reps" window the digest already flags as a possible misprint).
- **Ch. 43** (`prp-lightest-breather-universal`): `prp-an-linear-transmission` item 5 gives the classical factors that the coupling-independent lightest-breather factors must equal for all $a_r$ (cz2007 (7.5) for $a_2$); open problems: quantum meaning of the complex (charged) classes and of the real absorptions by charged defects; whether the missing charges can be reached by any defect process (classically they cannot, `prp-stored-charge` item 4).
