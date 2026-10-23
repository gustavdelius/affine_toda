# Referee prompt: half-line open problem 02
Tier: B
Brief: briefs/halfline_open_02_beyond_an.md
Report: briefs/halfline_open_02_report.md
Scripts: briefs/halfline_open_02_scripts/
- **Imports.** `dn_common.py` imports `code/boundary_common.py` read-only.
- **Core modules:** `dn_common.py`, `dn_scalar.py`, `proj_dn.py`, `dn_particle.py`, `dn_boot.py`, `caf.py`, `caf_all.py`, `fo_cmp.py`, `dn_classical.py`, `sym_vac2.py`, `spin_cl.py`, `dn_min.py`.
- **Small tests:** the `t*.py` files.
- **Superseded:** `sym_vac.py`, `t8.py`, `t9.py`, `fix*.py`.

Scope: step 1 ($d_n$, $\hat\epsilon=0$) and step 4 (the map) are done; step 3 ($c_2$) has only an analysis of its premises; step 2 ($d_n$, $s=\pm1$) was not started. Tier B because of a claimed error in a published paper (claim 7) and corrections to the brief (claim 10).

## Claims to referee

1. [Numerical] **Unitarity of the $\hat\epsilon=0$ matrix.** $K=A(\theta)J$ with $J$ antidiagonal, $\pm1$ entries, $\theta$-independent and $J^2=1$. Because $\bar V\cong V$ and $K$ is unique, unitarity is $A(\theta)A(-\theta)=1$.
2. [Numerical] **Crossing.**
   - The crossing matrix $C$ and the factor $c_u=\sin\pi\mu\sin\pi(\mu+T-\omega)/(\sin\pi(\mu+\omega)\sin\pi(\mu+T))$ are given.
   - The boundary-state vector lies in $L(2\omega_1)$, where $\check R=1$, so $\Phi\equiv1$.
   - The functional equations are therefore DG (3.8) with the $d_n$ factor $F$.
3. [Numerical] **The new vector–vector scalar factor $F$** (Gamma product, report §1). It is unitary, satisfies crossing $F(T-\mu)=F(\mu)c_u(-\mu)$ ($n=4,5,6$, ≤1e-7), and has the stated pole and zero structure.
4. [Numerical] **The boundary scalar factor $A$** (Gamma product). $A(0)=1$, and it has no $(-T/2)_I$ prefactor: adding one makes crossing fail by $-1$. Both equations hold for $n=4,5,6$ (≤1e-6).
5. [Numerical] **The breather.** The singlet projection gives exactly 1 times the singlet, and $K_{B_1}$ is DG (3.11) with $T=(n-1)\omega$ times four extra blocks (≤1e-6).
6. [Numerical] **The particle amplitude** $K_1=(h-\frac B2)(\frac h2-1)(\frac h2+1-\frac B2)/[(h-1)(1-\frac B2)(\frac h2)(\frac h2-\frac B2)]$, $h=2n-2$.
   - Unitarity, crossing-unitarity, all vector-like bulk boundary bootstraps and $1+2\to1$ hold for $n=4..7$ (exact block arithmetic).
   - $K_a$ equals Castro-Alvaredo–Fring's $\hat R_a$, (Rhat) with (am3).
7. [Literature correction] **CAF's spinor $\tilde R_s$ (their (dn2)) violates their own crossing relation (cuR) for odd $\ell$.**
   - For odd $n$, $\hat R_s$ satisfies crossing-unitarity only with $S_{ss'}$ in place of $S_{ss}$.
   - This is proposed as an Appendix D entry, and the odd-$n$ spinor amplitudes are left open.
8. [Numerical] **The classical identification: $\hat\epsilon=0$ is $C\equiv+1$ for $d_4$ and $d_5$.**
   - The diagram-automorphism-invariant vacua reduce to Tzitzéica ($d_4$, $e^{f(0)}=9/8$) and sinh-Gordon ($d_5$, $e^{f(0)}=4-2\sqrt2$). They are regular static kinks.
   - The linearized reflection about them equals $K_1(B\to0)$ (1e-8), and the spinors reflect diagonally.
   - That these are the vacua is not proved: a $d_4$ lattice minimization was poorly converged.
9. [Observation] **Fateev–Onofri.** $K_1$ equals FO (r1c) as a block identity at $w_0=w_r=x/2$, with $h\to2n-2$, $x\to B/2$. This is suggestive of the $b_r$ folding; the coupling map is unmatched.
10. [Correction to the brief] **Three premises of the brief are wrong.**
    - hep-th/9907072 and hep-th/0002213 are bulk reflection-amplitude papers, not boundary ones.
    - FO's $C_r$ family is the $C_{\rm short}=0$ branch, not the all-$\pm1$ branch.
    - Folding $a_3$ quantum amplitudes cannot give $c_2$ quantum amplitudes, because the $c_2$ S-matrix (DGZ) is not the restricted $a_3$ one.
11. [Literature map] **The step-4 table** (report §5). CAF's coverage is taken from their section list and verified only for $d_n$ and $a_2$.
12. [Text] **The proposed book text**: open problem 2, a new paragraph in `sec-dn-vector-k`, and a new `sec-dn-uniform`.

## Coordinator notes (suspected weak points)

- **Claim 7 is the most likely place for a convention error, so check it first.**
  - For odd $n$ the spinors $s$, $s'$ of $d_n$ are each other's *antiparticles*. Boundary crossing-unitarity for a non-self-conjugate particle relates $K_a$ to $K_{\bar a}$, and the bulk factor on the right may be $S_{a\bar a}=S_{ss'}$, not $S_{aa}$, depending on how the boundary state pairs particles.
  - Derive the correct diagonal crossing-unitarity for non-self-conjugate particles from Ghoshal–Zamolodchikov. Compare with the book's `eq-boundary-crossing-diagonal` and with DG's (3.17) for $a_n$, where particles are not self-conjugate either.
  - Read CAF's (cuR) and (dn2) in their TeX source. If CAF's own (cuR) uses the $a\bar a$ form, the "violation" may be the worker's misapplication. Settle this before anything goes into Appendix D.
- **Claim 4, the prefactor.** The book was once burned by a misquoted prefactor of the diagonal $a_n$ scalar factor, $-(\frac h4\omega)_I$. The worker drops DG's $(-T/2)_I$ prefactor for $d_n$, citing a sign from $g(\mu,0)$ in $F$, and normalizes $F(0)=-1$. Check the sign and normalization chain independently, including whether $F(0)=-1$ is the right normalization for the book's $d_n$ soliton S-matrix (Part VI).
- **Claim 8.**
  - Verify the reductions: that $\vec\phi=f(e_2-e_3)$ is consistent with the $d_4$ field equations and gives Tzitzéica; and similarly for $d_5$.
  - Verify that the boundary values satisfy the $C\equiv+1$ boundary condition in the book's normalization of the $d_n$ boundary potential, including the root-multiplicity factors $\sqrt{n_j}$, which differ between end and middle nodes for $d_n$.
  - Recompute the linearized $K$ about the $d_5$ kink yourself at one $\theta$.
  - Check whether a lower-energy configuration outside the invariant sector exists, using any minimizer you trust.
- **Claim 3.** Compare $F$ with the published $d_n$ vector soliton S-matrix scalar factor, for example Johnson 1997 or Gandenberger–MacKay–Watts, after converting conventions. Run the bulk breather test the worker did not do: the lowest breathers should reproduce $S_{11}$ of $d_n$.
- **Claim 6.** Recompute $K_1$ from the breather formula by the continuation $\omega\to-2/B$ yourself, and check unitarity and crossing-unitarity at one generic $\theta$ for $n=4$ and $n=5$.
- **Claim 10.** Check the three corrections against the arXiv abstracts and full text.
