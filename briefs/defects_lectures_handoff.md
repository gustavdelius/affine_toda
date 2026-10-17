# Handoff: the lecture notes "Windows in spacetime" (Lectures 6–9)

Written 2026-10-17 by the Claude session that drafted Lectures 1–5, for the Claude agent that writes Lectures 6–9. Read `AGENTS.md` first: its rules apply unchanged. This file adds what is specific to the notes, the state of the work, and a ready plan for Lecture 6. It can be given verbatim as the prompt of a fresh session started in this repository. Do not spawn subagents.

## 1. The task and the editor's decisions

On 2026-10-15 the editor (Gustav Delius) asked for lecture notes on Part VIII of the book (integrable defects). They are separate from the book and tell a compelling story: where each idea came from, and the decisive calculations written out in full, so that students can rediscover the results. Depth matters more than breadth. The editor approved this plan unchanged:

- a standalone Quarto book in `lecture-notes/defects/`, titled "Windows in spacetime", with 9 lectures of 12–18 pages each;
- sine-Gordon theory is the worked laboratory;
- $a_r^{(1)}$ appears in Lectures 5 and 9 only, but see §6.9 for the small $a_r$ section that Lecture 6 owes Lecture 5;
- the lecture list and the questions are in `lecture-notes/defects/index.qmd`.

| | Lecture | Question | State |
|---|---|---|---|
| 1 | Mirrors and windows | | merged (PR #30) |
| 2 | Bäcklund's machine | | merged |
| 3 | Freezing it at a point | | merged |
| 4 | What gets through | | merged |
| 5 | Why the cycle? | | merged |
| 6 | The classical algebra of a window | Why are the windows integrable, in the language of Lax pairs and r-matrices? | scripts and figures on branch `claude/defects-lecture-6`; text to write (§6) |
| 7 | Quantizing the window | Which symmetry survives in the quantum theory, and what does it fix? | not started (§7) |
| 8 | Windows, L-operators and Q-operators | What is a window wrapped around a cylinder? | not started (§7) |
| 9 | Folding, fusing and the classical breather | How are windows and mirrors related, and what is left of $\hbar$? | not started (§7) |

**Workflow, as the editor has used it.** The editor says "start Lecture N". Draft it, verify it, render it, and report. Commit or push only when asked. Each lecture went into a PR as one commit. The lecture notes PR #30 (Lectures 1–5) and the book PR #33 are both merged into `main`. The branch `claude/defects-lecture-6` holds the Lecture 6 scripts and this file.

**Improvements to the book.** Writing Lecture 5 produced improvements to book chapter 38. The editor asked for them in a separate PR from the lectures (#33, now merged). Do the same when a lecture finds something the book should say: offer it, and keep book changes out of the lecture PR.

## 2. Infrastructure of the notes (built in PR #30)

- **Book links.** `lecture-notes/book-links.lua` turns `[text](book:<id>)` into links to the published book. Never type a book number. Link only ids that exist on `main`: `python3 check-refs.py lecture-notes/defects` checks book links against the book's ids.
- **Check refs.** `python3 check-refs.py` (the book; it skips nested projects) and `python3 check-refs.py lecture-notes/defects` (the notes). Run both before a commit that touches `.qmd` files.
- **Render.** `./render-scratch.sh --project lecture-notes/defects` renders a scratch copy, in about 2–4 minutes, and lists unresolved references. Run it in the background. Never run `quarto render` in the checkout.
- **Scripts.**
  - They live in `code/lectures/`: `defects_lN_*.py` for the checks and `defects_lN_fM_*.py` for the figures. They are written from the lectures' own definitions, independently of the book's scripts.
  - Each prints one PASS or FAIL line per check, and exits with 1 on any FAIL.
  - List each script in `code/lectures/README.md`, in the columns script, lecture (ids), what it checks, and runtime.
  - Figure scripts import `code/figures/style.py` and write `lecture-notes/defects/figures/lN_*.svg`.
  - `script-links.lua` resolves script paths from the repository root. Cite a script at the claim, as a full path in inline code.
- **Bibliography.** Keys used only by the notes go in `lecture-notes/defects/references-notes.bib`, because the book's `references.bib` is nocite `@*`. Book keys can be cited directly.
- **Registering a lecture.** Add the file to `chapters:` in `lecture-notes/defects/_quarto.yml`. Update the status line in `index.qmd` ("Lectures 1 to 5 are written; ..."). Add the script rows to the README.
- **Visual check.**
  - Serve the scratch `_site` with `python3 -m http.server`.
  - Check wide displays at the 763 px text column (a 1280 px viewport). Split any equation wider than the column into an `aligned` block.
  - Kill the server by PID. Find the PID with `ss -ltnp | grep <port>`, because an old server from another session once survived its parent.

## 3. Voice and structure of a lecture

Re-read the previous lecture before drafting, for voice. Every lecture has this shape:

1. A title with `{#sec-...}`, then a callout `{.callout-note appearance="simple" icon=false title="In this lecture"}` with five or six bullets and a `*Prerequisites:*` line, which gives lecture sections and book links.
2. "The question" (`sec-lN-question`): the question in a blockquote, a `{.callout-tip title="Predict"}` box, and a `{.callout-note title="Where it comes from"}` box with the history. Every fact in the history box is checked against the source.
3. Sections that do the decisive calculations line by line, with `{.callout-note collapse="true" title="Check it yourself"}` boxes. Each box names the script and the item letters that check that section.
4. Two or three figures whose scripts check what they plot. Captions name the script.
5. "What we have learned" (`sec-lN-summary`), as bullets.
6. Exercises (`exr-lN-*`), each opening with *(Redo.)*, *(Rediscover.)* or *(Explore.)*. Hints go in `{.callout-tip collapse="true" title="Hint"}`.
7. "Where this comes from" (`sec-lN-sources`). It says who found what. It ends by listing what is "our presentation; we do not claim that it is new", and by linking the book's statements.

**Style.**

- Short declarative sentences. Define every symbol.
- Status labels as in the book: `**[Theorem]**` only with a complete proof in the text, `**[Established]**` for literature results the notes check. Every new [Theorem] has had no referee pass yet, so list it in the PR description.
- Never cite `briefs/` or `referee_scratch/`, and never mention a referee.
- At most one quotation per response from a source, under 15 words.
- Use they/them for people whose pronouns are not stated.

**Ids** are lecture-prefixed: `sec-l6-...`, `eq-l6-...`, `prp-l6-...`, `thm-l6-...`, `fig-l6-...`, `exr-l6-...`.

## 4. Conventions of the lectures

**Sine-Gordon (Lectures 2–4, 6–8).**

- The field $\Phi=\beta_{\rm SG}\phi$ is dimensionless. The action is $\frac1{\beta_{\rm SG}^2}\int[\frac12\Phi_t^2-\frac12\Phi_x^2-m_0^2(1-\cos\Phi)]$.
- $\partial_\pm=\partial_t\pm\partial_x$. The field equation is $\partial_+\partial_-\Phi=-m_0^2\sin\Phi$.
- The kink is $4\arctan e^{m_0x}$, with mass $M=8m_0/\beta_{\rm SG}^2$.
- The window has $\Phi$ on $x<0$ and $\Psi$ on $x>0$, with $\Sigma=\Phi+\Psi$ and $\Delta=\Phi-\Psi$ (not $A$, $B$).
- Type I$_-$ has $L_d=\Phi\Psi_t-\mathcal B$, with
  - $\mathcal B=-2m_0[\sigma(\cos\frac\Sigma2-1)+\sigma^{-1}(\cos\frac\Delta2-1)]$,
  - $P_d=2m_0[\sigma(\cos\frac\Sigma2-1)-\sigma^{-1}(\cos\frac\Delta2-1)]$,
  - $\sigma=e^{-\theta_d}$.
  - The light-cone pieces are $\mathcal B-P_d=-4m_0\sigma(\cos\frac\Sigma2-1)$ and $\mathcal B+P_d=-\frac{4m_0}\sigma(\cos\frac\Delta2-1)$.
- The defect conditions (`eq-l4-conditions`) are $\partial_+\Delta=-2\sigma m_0\sin\frac\Sigma2$ and $\partial_-\Sigma=\frac{2m_0}\sigma\sin\frac\Delta2$ at $x=0$. They are Lecture 2's Bäcklund transformation $\mathfrak B_a$ with $a=-\sigma$. Type I$_+$ is type I$_-$ followed by $\Psi\to-\Psi$.
- Solitons and the delay:
  - $Y=\exp(m_0(x\cosh\theta-t\sinh\theta)+c)$ and $\Psi=4\arctan(zY)$, with $z=-\coth\frac{\theta-\theta_d}2$;
  - linear waves cross with $\mathbb T(\theta)=z(\theta+\frac{i\pi}2)$;
  - the Bianchi coefficient is $\frac{a_1+a_2}{a_1-a_2}$ with $a=-e^{-\theta}$.
- A jump of $2\pi$ stores $M(\cosh\theta_d,\sinh\theta_d)$.

**Lecture 5 ($a_r^{(1)}$).**

- Real coupling. The fields $\vec\phi$, $\vec\psi$ are $\beta$ times the book's fields. $V=m^2\sum n_j(e^{\vec\alpha_j\cdot\vec\phi}-1)$.
- The kinetic matrices are sans-serif $\mathsf D,\mathsf E,\mathsf F,\mathsf P$, with $\mathsf D=1-\mathsf E$. The coordinates are $\vec\xi=\mathsf D^{\mathsf T}\vec\phi+\mathsf D\vec\psi$ and $\vec\eta=\vec\phi-\vec\psi$, and $\mathcal B=f(\vec\xi)+g(\vec\eta)$.
- In components, $\phi_k\in\mathbb R^h$ with $\sum_k\phi_k=0$:
  - $f_j=e^{\psi_j-\phi_{j+1}}$ and $g_j=e^{\phi_j-\psi_j}$;
  - the light-cone conditions are $\partial_+(\phi_k-\psi_k)=-m\sigma(f_k-f_{k-1})$ and $\partial_-(\psi_k-\phi_{k+1})=\frac m\sigma(g_k-g_{k+1})$;
  - $\omega=e^{2\pi i/h}$, and $\nu$ is used for frequency.

**New in Lecture 6. All of it is verified (§6).** The $2\times2$ matrices are sans-serif lowercase:

- $\mathsf h=\mathrm{diag}(1,-1)$, $\mathsf e=e_{12}$, $\mathsf f=e_{21}$, $\mathsf s=\mathsf e+\mathsf f$, and $\mathsf y=-i\mathsf e+i\mathsf f$. Do not use Pauli $\sigma$'s, which clash with the defect parameter.
- $\mathsf s_\Phi:=e^{i\Phi\mathsf h/4}\mathsf s\,e^{-i\Phi\mathsf h/4}=e^{i\Phi/2}\mathsf e+e^{-i\Phi/2}\mathsf f$. Then $\mathsf s_0=\mathsf s$, $\mathsf s_{2\pi}=-\mathsf s$, $\partial_\Phi\mathsf s_\Phi=\frac i4[\mathsf h,\mathsf s_\Phi]$ and $[\mathsf s_\Phi,\mathsf s_{-\Phi}]=2i\sin\Phi\,\mathsf h$.
- The Hadamard matrix $\mathsf u=\frac1{\sqrt2}\begin{pmatrix}1&1\\1&-1\end{pmatrix}$ satisfies $\mathsf{usu}=\mathsf h$, $\mathsf{uhu}=\mathsf s$ and $\mathsf{uyu}=-\mathsf y$.
- Use the book's names $\mathcal A_\pm$, $a_t$, $a_x$, $T$, $\mathcal T_d$ and $\tau_d$.

## 5. Process lessons from Lectures 1–5

- Rule 1 bites. Every lecture found slips in its own first drafts: log branches in soliton formulas, a false conjugation identity, a wrong attribution. Write the script before polishing the text. Mutation-test it: break a sign or a factor, and the expected check must FAIL.
- **Sympy pitfalls.**
  - Use `.rewrite(sp.exp)` before `simplify` for trigonometric identities.
  - Never call `.rewrite` on a Python int.
  - `sp.solve` on exact expressions with `exp(I*rational)` can silently find no solution. Use a numeric least-squares fit at a random point instead.
- Regularized transport matrices lose precision when products are formed in a non-diagonal basis. Rotate first, and take only the dominant entry.
- The session hook forbids writing to other worktrees' files. To build a commit for another branch, use a temporary `GIT_INDEX_FILE` with `read-tree`, `apply --cached`, `write-tree` and `commit-tree`.
- Literature: every attribution in the notes was checked against the arXiv abstract or the paper. Keep doing it. The book's chapters cite equation numbers of the sources, which were checked during the book's referee passes, and the notes may reuse them.

## 6. Lecture 6, "The classical algebra of a window": the plan, ready to write

File: `lecture-notes/defects/06-classical-algebra.qmd`. The scripts are already on this branch:

- `code/lectures/defects_l6_lax.py`: 16 checks, (a)–(o), about 80 s;
- `defects_l6_f1_holonomy.py`, which writes `l6_holonomy.svg`;
- `defects_l6_f2_charges.py`, which writes `l6_charges.svg`.

The docstring of `defects_l6_lax.py` states every check. Add the README rows with the final ids, and the lecture to `_quarto.yml` and the `index.qmd` status line.

**Storyline.** Lecture 3 found the window from momentum, which is one charge. Integrability needs infinitely many charges, in involution. The tool is the Lax pair. The window is found a third time, as a frozen gauge transformation: first by Bäcklund in Lecture 2, then by momentum in Lecture 3, now by the Lax pair. Its matrix is a classical L-operator for exactly the bulk r-matrix, and the kinetic term that momentum fixed is exactly the Poisson bracket the r-matrix requires.

**Predict box.** A boundary has a field-independent matrix $K$. Should the window's matrix depend on the fields? On $\lambda$? What should happen to it at the defect rapidity?

**Where it comes from.** Each item below was checked against its arXiv abstract on 2026-10-17, or for the last two against the bibliographic record:

- Bowcock, Corrigan and Zambon built the Lax pair of the defect from two Lax pairs on overlapping regions [@bcz2004, Sec. 4; @bcz2004b, Sec. 3], following the book's attribution. hep-th/0401020 constructs a Lax pair and finds that only $a_r^{(1)}$ allows integrable defects.
- Caudrelier (0704.2326, IJGMMP 5 (2008) 1085) used inverse scattering. He built the generating function of all the modified charges, identified the defect's contribution to all orders through a defect matrix, and observed that the defects are Bäcklund transformations localized at a point. He also classified the defect matrices, and noted that the r-matrix approach needs a lattice regularization.
- Habibullin and Kundu (0709.4611) used a lattice regularization. They found the R-matrices and q-algebra, and all charges algorithmically.
- Avan and Doikou (1205.1661 and 1110.4728) took an algebraic route. The r-matrix algebra comes first, and "sewing conditions" follow from continuity of the time components.
- Caudrelier (1411.5171) and Caudrelier–Kundu (1411.0418) used a second, multisymplectic Poisson bracket. In it the defect conditions are a canonical transformation, and the problem of fields at coinciding points disappears.
- Sklyanin 1988: the double row at a boundary (`sklyanin1988`).
- For the bulk, cite these three, to be added to `references-notes.bib` (§6.10):
  - Lax's pairs of 1968;
  - Ablowitz–Kaup–Newell–Segur 1973, who solved sine-Gordon by inverse scattering;
  - Faddeev–Takhtajan's book, for the Hamiltonian theory and the classical r-matrix.

### 6.1 The Lax pair (`sec-l6-lax`)

$$
\mathcal A_+=\tfrac i4\partial_+\Phi\,\mathsf h+\tfrac{m_0\lambda}2\mathsf s_\Phi,\qquad
\mathcal A_-=-\tfrac i4\partial_-\Phi\,\mathsf h-\tfrac{m_0}{2\lambda}\mathsf s_{-\Phi}.
$$

- **Zero curvature, worked by hand.** $\partial_+\mathcal A_--\partial_-\mathcal A_++[\mathcal A_+,\mathcal A_-]=-\frac i2(\partial_+\partial_-\Phi+m_0^2\sin\Phi)\mathsf h$.
  - The $\lambda$ terms cancel because $\partial_-\mathsf s_\Phi=\frac i4\partial_-\Phi[\mathsf h,\mathsf s_\Phi]$. Likewise for $1/\lambda$.
  - $\lambda^0$ gives $-\frac i2\partial_+\partial_-\Phi\,\mathsf h-\frac{m_0^2}4[\mathsf s_\Phi,\mathsf s_{-\Phi}]$.
  - [checked: (a)]
- **The book's pair.** This is the book's `eq-lax-pair` for $a_1^{(1)}$ at imaginary coupling, with $\sqrt2\beta_b\phi=\Phi$, $m=m_0/2$ and $H=\mathsf h/\sqrt2$. [(a)]
- **Boosts.** A boost by $\vartheta$ is $\lambda\to e^\vartheta\lambda$. [(a)]
- **The two components.**
  - $a_t=\frac12(\mathcal A_++\mathcal A_-)=\frac i4\Phi_x\mathsf h+\frac{m_0}4(\lambda\mathsf s_\Phi-\lambda^{-1}\mathsf s_{-\Phi})$;
  - $a_x=\frac i4\Phi_t\mathsf h+\frac{m_0}4(\lambda\mathsf s_\Phi+\lambda^{-1}\mathsf s_{-\Phi})$;
  - the linear problem is $(\partial_x+a_x)\chi=(\partial_t+a_t)\chi=0$.
- **The vacuum.**
  - $a_x(0)=\frac{m_0}4(\lambda+\lambda^{-1})\mathsf s$ and $a_t(0)=\frac{m_0}4(\lambda-\lambda^{-1})\mathsf s$.
  - At $\lambda=-e^{-\theta}$, $\chi=c_+e^{m_0\xi/2}(1,1)+c_-e^{-m_0\xi/2}(1,-1)$ with $\xi=x\cosh\theta-t\sinh\theta$.
  - The ratio of the two components is $\frac{c_+}{c_-}e^{m_0\xi}$, the soliton's $Y$. So $\lambda=-e^{-\theta}$ is a rapidity, as Lecture 2's $a=-e^{-\theta}$ was. [(b)]
- **Transport.** $\partial_xT(x,y)=-a_x(x)T$, and flatness gives $\partial_tT=-a_t(x)T+Ta_t(y)$. Two lines: $K=\partial_tT+a_t(x)T-Ta_t(y)$ obeys $\partial_xK=-a_xK$ with $K(y,y)=0$.

### 6.2 Why Bäcklund's machine works (`sec-l6-gauge`), a [Theorem], proof in full

Write $D=e^{i\Delta\mathsf h/4}$ and $O=e^{-i\Sigma\mathsf h/4}\mathsf s$, and set $L=D+\frac\sigma\lambda O$.

**Claim.** $\partial_\pm L=L\mathcal A_\pm(\Psi)-\mathcal A_\pm(\Phi)L$ if and only if $\Psi=\mathfrak B_{-\sigma}\Phi$, that is $\partial_+\Delta=-2\sigma m_0\sin\frac\Sigma2$ and $\partial_-\Sigma=\frac{2m_0}\sigma\sin\frac\Delta2$. [(c)]

So $\mathcal A(\Psi)=L^{-1}\mathcal A(\Phi)L+L^{-1}\partial L$, and the curvature transforms by conjugation. This explains the "surprise" of Lecture 2, that the Bäcklund system is consistent.

**Proof of the $\partial_+$ equation, order by order.**

- $\lambda^1$: $D\mathsf s_\Psi=\mathsf s_\Phi D$, an identity.
- $\lambda^0$: $\partial_+D=\frac i4\partial_+\Delta\,\mathsf hD$, and $O\mathsf s_\Psi-\mathsf s_\Phi O=-2i\sin\frac\Sigma2\,\mathsf hD$, give exactly the first equation.
- $\lambda^{-1}$: an identity, because $O\mathsf h=-\mathsf hO$.

**The $\partial_-$ equation.**

- $\lambda^0$: an identity.
- $\lambda^{-1}$: $D\mathsf s_{-\Psi}-\mathsf s_{-\Phi}D=2i\sin\frac\Delta2\,\mathsf hO$ gives the second equation.
- $\lambda^{-2}$: $O\mathsf s_{-\Psi}=\mathsf s_{-\Phi}O$.

$\det L=1-\sigma^2/\lambda^2$, and $L$ depends on $\sigma/\lambda$ only, which is $-e^{\theta-\theta_d}$ at $\lambda=-e^{-\theta}$. [(e), (n)]

### 6.3 Freezing it (`sec-l6-frozen`), a [Theorem] and the decisive calculation

At $x=0$ only $\partial_t$ is meaningful. Gluing $\chi_\Phi(0)=G\chi_\Psi(0)$ is compatible with the time evolution if and only if $\partial_tG=Ga_t(\Psi)-a_t(\Phi)G$.

**Ansatz.** $G=G_0+\lambda^{-1}G_1$.

- Order $\lambda^1$: $G_0\mathsf s_\Psi=\mathsf s_\Phi G_0$, which forces $G_0=pD+q\bar O$, where $\bar O:=\mathsf s\,e^{-i\Sigma\mathsf h/4}=e^{i\Sigma\mathsf h/4}\mathsf s$.
- Order $\lambda^{-2}$: $G_1\mathsf s_{-\Psi}=\mathsf s_{-\Phi}G_1$, which forces $G_1=p'D^{-1}+q'O$.
- These are the general solutions of the two linear equations. [(d)]

**The block identities.**

- The extreme orders are identities: $\bar O\mathsf s_\Psi=\mathsf s_\Phi\bar O$ and $D^{-1}\mathsf s_{-\Psi}=\mathsf s_{-\Phi}D^{-1}$.
- The middle orders use:
  - $O\mathsf s_\Psi-\mathsf s_\Phi O=-2i\sin\frac\Sigma2\,\mathsf hD$;
  - $D^{-1}\mathsf s_\Psi-\mathsf s_\Phi D^{-1}=-2i\sin\frac\Delta2\,\mathsf h\bar O$;
  - $D\mathsf s_{-\Psi}-\mathsf s_{-\Phi}D=2i\sin\frac\Delta2\,\mathsf hO$;
  - $\bar O\mathsf s_{-\Psi}-\mathsf s_{-\Phi}\bar O=2i\sin\frac\Sigma2\,\mathsf hD^{-1}$.
- These were derived by hand in the session. Item (d) checks the resulting four equations. Consider adding the eight block identities as explicit checks.

**The middle orders.** Separate diagonal from off-diagonal, and the identity from $\mathsf h$. The orders $\lambda^0$ and $\lambda^{-1}$ are equivalent to $p,q,p',q'$ constant together with

1. $p\,\partial_+\Delta=-2m_0q'\sin\frac\Sigma2$,
2. $p'\,\partial_-\Delta=2m_0q\sin\frac\Sigma2$,
3. $q\,\partial_+\Sigma=-2m_0p'\sin\frac\Delta2$,
4. $q'\,\partial_-\Sigma=2m_0p\sin\frac\Delta2$.

Item (d) checks that the 8 middle coefficients are an invertible combination of these four equations and the four time derivatives.

**The pairing argument, which needs no case enumeration.** The defect values and the velocities are free at one instant (Lecture 3).

- (1) and (4) pair $p$ with $q'$: if one of them is zero and the other is not, one equation reads "$0=$ constant times a sine", a relation among the values.
- (2) and (3) pair $q$ with $p'$ in the same way.
- If all four constants are nonzero, the four equations fix all four first derivatives, so the velocities are not free. Item (d) checks the rank 4.

So either $q=p'=0$, which is type I$_-$ with $\sigma=q'/p$ and $G=p\,L_d$, or $p=q'=0$, which is type I$_+$ with $\sigma=p'/q$ and $G=q\,L_d^+$. Here $L_d^+=\mathsf s\,e^{-i\Sigma\mathsf h/4}+\frac\sigma\lambda e^{-i\Delta\mathsf h/4}=L_d(\Phi,-\Psi)\,\mathsf s$, and its conditions are those of I$_-$ with $\Psi\to-\Psi$. Underlying this, $\mathcal A_\pm(-\Psi)=\mathsf s\mathcal A_\pm(\Psi)\mathsf s$; consider adding it as a check. [(d), (n)]

**The result.** The window is the frozen gauge transformation:

$$
L_d(\lambda)=e^{i\Delta\mathsf h/4}+\tfrac\sigma\lambda e^{-i\Sigma\mathsf h/4}\mathsf s=\begin{pmatrix}e^{i\Delta/4}&\frac\sigma\lambda e^{-i\Sigma/4}\\\frac\sigma\lambda e^{i\Sigma/4}&e^{-i\Delta/4}\end{pmatrix}.
$$

- The $\Delta$-half of $\mathcal B$ is on the diagonal. The $\Sigma$-half, with $\sigma$, is off the diagonal.
- It equals $e^{i\Delta\mathsf h/4}(1+\frac\sigma\lambda\mathsf s_{-\Phi})$, the book's `eq-defect-lax-matrix` form; consider adding this as an explicit check.
- At $\lambda=-\sigma$ it has rank 1, with kernel $(1,e^{i\Phi/2})$. [(e)]
- Remark: $L_d(\sigma)(\Phi,\Psi)\,L_d(-\sigma)(\Psi,\Phi)=\det L_d$, which is $\mathfrak B_{-a}=\mathfrak B_a^{-1}$ of Lecture 2. [(n)]

### 6.4 What the matrix knows (`sec-l6-knows`)

- On the vacuum, $L_d=1+\frac\sigma\lambda\mathsf s$, with eigenvalues $1\pm\sigma/\lambda$ on $(1,\pm1)$.
- At $\lambda=-e^{-\theta}$, the ratio $(1-\sigma/\lambda)/(1+\sigma/\lambda)$ is Lecture 4's delay $z=-\coth\frac{\theta-\theta_d}2$. At $\lambda=ie^{-\theta}$ it is the transmission factor $\frac{\sinh\theta_d+i\cosh\theta}{\cosh\theta_d-i\sinh\theta}$, with modulus 1. [(e)]
- Heuristic reason, to be phrased as such: $Y$ is the ratio of the vacuum components along $(1,1)$ and $(1,-1)$, and $\chi_\Psi=L_d^{-1}\chi_\Phi$ multiplies them by the inverse eigenvalues, so $Y\to zY$. The exact argument goes through the soliton's bound state at $\lambda=-e^{-\theta}$ and its norming constant. Use these facts: in the original gauge, $a_x(2\pi)=-a_x(0)$, so the decaying solution is along $(1,1)$ at both ends; $L_d(0,0)=1+\frac\sigma\lambda\mathsf s$ and $L_d(2\pi,2\pi)=1-\frac\sigma\lambda\mathsf s$ [(n)]. Leave the details as an Explore exercise. They are not fully checked, so do not state them as a theorem.
- $\det L_d=0$ at $\lambda=\pm\sigma$, that is at $\theta=\theta_d$ (absorption) and at $\theta_d+i\pi$.

### 6.5 The generating function (`sec-l6-monodromy`), a [Theorem]

**On a circle.** Put $\Phi$ on $[-\ell,0)$ and $\Psi$ on $(0,\ell]$, with $\pm\ell$ identified, and the fields continuous there and without winding.

- $\mathcal T_d=T^\Psi(\ell,0)L_d^{-1}T^\Phi(0,-\ell)$.
- $\partial_t\mathcal T_d=-a_t(\Psi;\ell)\mathcal T_d+\mathcal T_da_t(\Phi;-\ell)=[\mathcal T_d,a_t(\ell)]$, so $\tau_d=\operatorname{tr}\mathcal T_d$ is conserved.
- The time equation of $L_d^{-1}$ is $\partial_tL_d^{-1}=L_d^{-1}a_t(\Phi)-a_t(\Psi)L_d^{-1}$.

**On the line.** Use the gauge $\chi=e^{i\Phi\mathsf h/4}\tilde\chi$, in which every vacuum $2\pi n$ has the same $\tilde a$, namely $\frac{m_0}4(\lambda\pm\lambda^{-1})\mathsf s$, and $\tilde L_d=1+\frac\sigma\lambda\mathsf s_{-\Sigma}$ [(e)]. The regularized monodromy $e^{k\mathsf s x_+}\tilde T_+\tilde L_d^{-1}\tilde T_-e^{-k\mathsf s x_-}$, with $k=\frac{m_0}4(\lambda+\lambda^{-1})$, has constant diagonal entries in the eigenbasis of $\mathsf s$. [(f), numerically along an exact passage; fails with $L_d\to1$ or $1.3\sigma$]

**Figure `fig-l6-holonomy`** (`l6_holonomy.svg`, by `defects_l6_f1_holonomy.py`, 3 checks).

- (a) The passage with $\theta=0.6$, $\theta_d=1$, shown by its charge density, and a rectangle crossing the window, with $L_d^{-1}$ at $t_1$ and $L_d$ at $t_2$.
- (b) $\Vert H-1\Vert$ against $t_2$, for $\lambda=0.8+0.5i$. With $L_d$ it is about $10^{-11}$; with $L_d\to1$ or $\sigma\to1.3\sigma$ it is of order 1 to 30 once the soliton arrives.
- The caption must say that flatness with a seam means a holonomy of 1.

**Contrast with the mirror.** Sklyanin's double row (book `eq-double-row`) contains $T(\lambda)$ and $T(1/\lambda)$, so only $Q_s+Q_{-s}$ survive. The window has one $\lambda$, so $Q_s$ and $Q_{-s}$ survive separately, closing Lecture 1's remark in `sec-l1-momentum`.

### 6.6 The charges (`sec-l6-charges`), the decisive calculation

Present it on the circle. On the line the window terms are the same, checked in (g) and (h) on the line.

**Steps at $\lambda\to\infty$.**

1. In the gauge $e^{i\Phi\mathsf h/4}$, $\tilde a_x=\frac i4\partial_+\Phi\,\mathsf h+\frac{m_0\lambda}4\mathsf s+\frac{m_0}{4\lambda}\mathsf s_{-2\Phi}$.
2. Rotating by $\mathsf u$ gives $\hat a=\frac{m_0\lambda}4\mathsf h+\frac i4\partial_+\Phi\,\mathsf s+\frac{m_0}{4\lambda}(\cos\Phi\,\mathsf h-\sin\Phi\,\mathsf y)$ and $\hat L_d=1+\frac\sigma\lambda(\cos\frac\Sigma2\mathsf h-\sin\frac\Sigma2\mathsf y)$.
3. Write $\hat T(x,y)=(1+W(x))e^{Z(x,y)}(1+W(y))^{-1}$, with $W$ off-diagonal and $Z$ diagonal. Then
   - $Z'=-\mathrm{diag}\,\hat a-(\mathrm{off}\,\hat a\,W)_{\rm diag}$;
   - $W'=-\mathrm{off}\,\hat a-[\mathrm{diag}\,\hat a,W]+W\,\mathrm{off}\,\hat a\,W$.
4. $W_1=\frac{\partial_+\Phi}{2m_0}\mathsf y$, using $[\mathsf h,\mathsf y]=-2i\mathsf s$ and $\mathsf{sy}=i\mathsf h$.
5. $Z'_{22}=\frac{m_0\lambda}4+\frac1\lambda\big[\frac{m_0}4-\frac1{4m_0}(\frac12(\partial_+\Phi)^2+V)\big]+O(\lambda^{-2})$. The 22 entry dominates for large positive $\lambda$.
6. $\ln\tau_d=Z^\Phi_{22}+Z^\Psi_{22}+\ln N_{22}$, with $N=(1+W^\Psi(0))^{-1}\hat L_d^{-1}(1+W^\Phi(0))$, and $\ln N_{22}=\frac\sigma\lambda\cos\frac\Sigma2+O(\lambda^{-2})$.
7. The $\lambda^{-1}$ coefficient is $\frac{m_0\ell}2+\sigma-\frac1{4m_0}[\beta^2(E-P)+\mathcal B-P_d]$, so Lecture 3's modified momentum reappears. [(g)]

**At $\lambda\to0$.**

- In the gauge $e^{-i\Phi\mathsf h/4}$: $a'_x=\frac i4\partial_-\Phi\,\mathsf h+\frac{m_0}{4\lambda}\mathsf s+\frac{m_0\lambda}4\mathsf s_{2\Phi}$, and $L'=e^{i\Delta\mathsf h/2}+\frac\sigma\lambda\mathsf s$.
- In the rotated basis $\hat L'=\frac\sigma\lambda\mathsf h+\cos\frac\Delta2+i\sin\frac\Delta2\,\mathsf s$, and $\ln$ of its dominant inverse entry is $\ln(-\lambda/\sigma)+\frac\lambda\sigma\cos\frac\Delta2+O(\lambda^2)$.
- So $E+P+(\mathcal B+P_d)/\beta^2$ is conserved. [(h)]
- Here $W'_1=\frac{\partial_-\Phi}{2m_0}\mathsf y$ and $W'_2=-\frac i{2m_0^2}\partial_-^2\Phi\,\mathsf s$, whereas $W_2=+\frac i{2m_0^2}\partial_+^2\Phi\,\mathsf s$.

**Spin $\mp3$** (all exact, symbolic and numerical, in (g) and (h)).

- The bulk density at $\lambda^{-3}$ is $\frac1{8m_0^3}q_-$ with $q_-=(\partial_+^2\Phi)^2-\frac14(\partial_+\Phi)^4+m_0^2\cos\Phi(\partial_+\Phi)^2$, minus $\partial_x(\partial_+\Phi\,\partial_+^2\Phi)/(4m_0^3)$. This is the book's `eq-spin3-current` continued to sine-Gordon.
- The window term, after the boundary terms and the first defect condition, which makes the $\partial^2$ terms cancel:
$$
\mathcal D_{-3}=\frac\sigma{2m_0^2}\cos\tfrac\Sigma2\,\partial_+\Phi\,\partial_+\Psi+\sigma^3\big(\cos\tfrac\Sigma2-\tfrac23\cos^3\tfrac\Sigma2\big).
$$
- At $\lambda\to0$, $\mathcal J'=+\partial_-\Phi\,\partial_-^2\Phi/(4m_0^3)$ and
$$
\mathcal D_{+3}=-\frac1{2m_0^2\sigma}\cos\tfrac\Delta2\,\partial_-\Phi\,\partial_-\Psi+\sigma^{-3}\big(\cos\tfrac\Delta2-\tfrac23\cos^3\tfrac\Delta2\big).
$$
  Note the minus sign.
- The conserved charges are $\frac1{8m_0^3}[\int_{x<0}q_\mp(\Phi)+\int_{x>0}q_\mp(\Psi)]+\mathcal D_{\mp3}$.
- A soliton carries $\beta^2(E\mp P)=8m_0e^{\mp\theta}$, and $-\frac23e^{\mp3\theta}$ for spin $\mp3$. The odd jump $(0,2\pi)$ stores exactly these values at $\theta_d$; subtract the even-vacuum constant $\sigma^{\pm3}/3$. [`defects_l6_f2_charges.py`; add an explicit check of $-\frac23e^{\mp3\theta}$ for general $\theta$]
- Do not make a labelled [Theorem] of the expansion. It is a formal asymptotic expansion, as in the book's `prp-spin-two`(2), which is labelled "sketch". Present the window terms as the result of a calculation checked by the script.

**Figure `fig-l6-charges`** (`l6_charges.svg`, 5 checks).

- A passage with $\theta_d=0.8$ and $\theta=\theta_d-0.01$, which is held for about 6 time units.
- Four panels, for spins $-1,+1,-3,+3$. Each shows the share of the field on $x<0$, of the window and of the field on $x>0$, and their sum.
- The window holds the whole charge during the hold.
- The spin $\pm3$ shares are not positive: they range over about $-0.2$ to $1.2$, because the densities are not positive, and only the sum is fixed. Say so in the caption.

### 6.7 The Poisson bracket of a window (`sec-l6-poisson`), a [Theorem] with a hand proof

- **The defect bracket.** $L_d/\beta^2=(\Phi\Psi_t-\mathcal B)/\beta^2$ is $p\dot q-H$ with $q=\Psi$ and $p=\Phi/\beta^2$, so $\{\Phi,\Psi\}=-\beta^2$ for type I$_-$ and $+\beta^2$ for type I$_+$. The convention $\{\Phi(x),\Phi_t(y)\}=\beta^2\delta$ matches the book's `eq-dc-defect-bracket` with $\Phi=\beta\phi$. [(j)]
- **The bulk r-matrix**, [Established], with the book's `prp-sklyanin` and Faddeev–Takhtajan:
  - $\{a_x(\lambda_1;x)\overset\otimes,a_x(\lambda_2;y)\}=\delta(x-y)[r(\lambda_1/\lambda_2),a_x(\lambda_1)\otimes1+1\otimes a_x(\lambda_2)]$;
  - $r(z)=-\frac{\beta^2}{16(z^2-1)}[(z^2+1)\mathsf h\otimes\mathsf h+4z\Pi]$ with $\Pi=\mathsf e\otimes\mathsf f+\mathsf f\otimes\mathsf e$;
  - $\{T_1,T_2\}=-[r,T_1T_2]$.
  - It is the book's $r$ for $h=2$ with $b^2=-\beta^2/2$. It satisfies the classical Yang–Baxter equation, $r^{t_1t_2}(z)=-r(1/z)$, and the lattice sign check $c\to-1$. [(i)] Its derivation is a Redo exercise.
- **Hand proof of $\{L_1,L_2\}=[r,L_1L_2]$.**
  - $\{\Delta,\Sigma\}=2\{\Phi,\Psi\}$. With $\kappa=\{\Phi,\Psi\}/8$: $\{D_1,O_2\}=\kappa\,\mathsf h\otimes\mathsf h\,D_1O_2$, $\{O_1,D_2\}=-\kappa\,\mathsf h\otimes\mathsf h\,O_1D_2$, and the brackets $DD$ and $OO$ vanish. Hence $\{L_1,L_2\}=\kappa\sigma\,\mathsf h\otimes\mathsf h(\lambda_2^{-1}D_1O_2-\lambda_1^{-1}O_1D_2)$.
  - $r$ commutes with $D_1D_2$ and with $O_1O_2$.
  - $[\mathsf h\otimes\mathsf h,X]=2\,\mathsf h\otimes\mathsf h\,X$ for $X=D_1O_2$ and $O_1D_2$.
  - $[\Pi,D_1O_2]=-\mathsf h\otimes\mathsf h\,O_1D_2$ and $[\Pi,O_1D_2]=-\mathsf h\otimes\mathsf h\,D_1O_2$. Verify by hand with $\mathsf{ef}=\frac12(1+\mathsf h)$ and $\mathsf{fe}=\frac12(1-\mathsf h)$. [(j)]
  - With $r=\alpha\,\mathsf h\otimes\mathsf h+\gamma\Pi$ and $z=\lambda_1/\lambda_2$, matching gives $2\alpha-\gamma/z=\kappa$ and $2\alpha-\gamma z=-\kappa$. So $r=\frac\kappa{2(z^2-1)}[(z^2+1)\mathsf h\otimes\mathsf h+4z\Pi]$: **the window by itself determines the trigonometric r-matrix**, and it is the bulk one exactly when $\kappa=-\beta^2/8$, that is $\{\Phi,\Psi\}=-\beta^2$, the bracket of the kinetic term that momentum fixed. Type I$_+$ needs $+\beta^2$, which its kinetic term gives. [(j)]
- **Involution, a corollary.**
  - $\{L^{-1}_1,L^{-1}_2\}=-[r,L^{-1}_1L^{-1}_2]$, the same algebra as $T$.
  - Products of Poisson-commuting factors with this algebra have it too, so $\mathcal T_d$ does, and $\{\tau_d(\lambda_1),\tau_d(\lambda_2)\}=0$. [(k)]
  - A chain of windows with independent pairs has commuting traces, which makes it a classical integrable lattice; this points to Lecture 8. Mixing $L$ and $L^{-1}$ fails. [(k)]
- **The caveat.** State honestly, following the book's `sec-defect-r-matrix`, that the involution assumes the defect values Poisson-commute with the bulk. The equal-time bracket does not restrict to the boundary values of the fields (Caudrelier 2015, Sec. 1). The book's paragraph "What the proposition does not settle" has the literature.

### 6.8 Space as time and two windows (short sections)

**Space as time** (`sec-l6-space-time`). For $G=\int dt(\Phi\Psi_t-\mathcal B)$:

- $\delta G/\delta\Phi(0,t)=\Psi_t-\mathcal B_\Phi=\Phi_x$ and $-\delta G/\delta\Psi(0,t)=\Phi_t+\mathcal B_\Psi=\Psi_x$;
- so $-G$ generates the canonical map $(\Phi,-\Phi_x)\to(\Psi,-\Psi_x)$ in the picture where $x$ is time [(m)];
- this is Caudrelier 2015, Prop. 3.1, via the book. In that bracket there is no coinciding-point problem, and Caudrelier proved Liouville integrability [Caudrelier 2015, Sec. 3.3, as cited by the book].

**Two windows** (`sec-l6-two-windows`).

- $L(\sigma_1;\Phi_0,\Phi_1)L(\sigma_2;\Phi_1,\Phi_{12})=L(\sigma_2;\Phi_0,\Phi_2)L(\sigma_1;\Phi_2,\Phi_{12})$ holds for all $\lambda$ when $\Phi_{12}$ is given by Bianchi's formula `eq-l2-bianchi` with $a_k=-\sigma_k$, and a perturbed $\Phi_{12}$ fails [(l)]. Both sides flip sign together under a change of branch.
- This is Bianchi's square as a refactorization of gauge matrices: two windows can be exchanged, with the field between them changed. It is the classical counterpart of the window–window relations of Lecture 7.
- Do not claim novelty. Refactorization of Darboux matrices is classical, but its origin was not traced.

### 6.9 The $a_r^{(1)}$ window (`sec-l6-ar`), owed to Lecture 5

Lecture 5 says "Lecture 6 constructs the Lax pair of the windows of this lecture", so give a short section.

- In Lecture 5's real-coupling components, use the defining representation: $H\cdot v=\mathrm{diag}(v_k)$, $E_{\vec\alpha_j}=e_{j,j+1}$ and $E_{-\vec\alpha_j}=e_{j+1,j}$, with indices mod $h$.
- The Lax pair is $\mathcal A_+=\frac12\partial_+\vec\phi\cdot H+m\lambda\sum_je^{\vec\alpha_j\cdot\vec\phi/2}E_{\vec\alpha_j}$ and $\mathcal A_-=-\frac12\partial_-\vec\phi\cdot H-\frac m\lambda\sum_je^{\vec\alpha_j\cdot\vec\phi/2}E_{-\vec\alpha_j}$. Its curvature is minus the Toda equation `eq-l5-toda` on the diagonal.
- The defect matrix is $(L_d)_{kk}=g_k^{1/2}$ and $(L_d)_{j+1,j}=\frac\sigma\lambda f_j^{1/2}$, with all other entries zero. **The diagonal carries the $g$-half and the subdiagonal the $f$-half.**
- It satisfies the time equation on Lecture 5's light-cone conditions, for $h=2,\dots,5$ and complex $\lambda$, and fails for $1.2\sigma$. [(o)]
- It is the book's $e^{H\cdot(\phi-\psi)/2}[1+\frac\sigma\lambda\sum_je^{\vec\alpha_j\cdot\vec w/2}E_{-\vec\alpha_j}]$, because $\vec\alpha_j\cdot\vec w=\psi_j+\psi_{j+1}-2\phi_{j+1}$. That identity was derived by hand and is not yet in the script.
- Link the book's `prp-defect-lax`, `prp-defect-involution` and `prp-spin-two`: both spin $\pm2$ charges survive, one at each end of the $\lambda$-plane.

### 6.10 References, exercises, labels

**References to add to `references-notes.bib`:**

- `lax1968`: P. D. Lax, "Integrals of nonlinear equations of evolution and solitary waves", Comm. Pure Appl. Math. 21 (1968) 467–490.
- `akns1973`: M. J. Ablowitz, D. J. Kaup, A. C. Newell and H. Segur, "Method for solving the sine-Gordon equation", Phys. Rev. Lett. 30 (1973) 1262–1264, doi 10.1103/PhysRevLett.30.1262.
- `faddeev_takhtajan1987`: L. D. Faddeev and L. A. Takhtajan, *Hamiltonian Methods in the Theory of Solitons*, Springer, Berlin, 1987.
- The Caudrelier, Habibullin–Kundu, Avan–Doikou and BCZ keys are already in the book's `references.bib`.

**Exercises.**

- `exr-l6-flat` (Redo): zero curvature, the time evolution of $T$, boosts.
- `exr-l6-backlund` (Redo): the $\partial_-$ half of 6.2, and the gauge covariance of the curvature.
- `exr-l6-types` (Rediscover): the type-I$_+$ branch and $L_d^+=L_d(\Phi,-\Psi)\mathsf s$. On the vacuum it is $\mathsf s+\sigma/\lambda\propto1+\frac\lambda\sigma\mathsf s$, the book's $\tilde L_{\rm I_+}$ for $h=2$.
- `exr-l6-zero` (Redo): the $\lambda\to0$ expansion.
- `exr-l6-norming` (Explore): the delay from the bound state.
- `exr-l6-rmatrix` (Redo): derive $r$, and the sign of $\{T,T\}$.
- `exr-l6-identities` (Redo): the four commutator identities and the solution for $\alpha,\gamma$.
- `exr-l6-chain` (Explore): a chain of windows as a lattice.
- `exr-l6-spin3` (Explore): derive $\mathcal D_{-3}$ and show directly that the flux of $q_-$ through $x=0$ is $-\frac d{dt}\mathcal D_{-3}$. Bowcock, Corrigan and Zambon found the window from the survival of this charge [@bcz2004]. Do not claim that they wrote $\mathcal D_{-3}$ out.
- `exr-l6-two` (Explore): two windows and Bianchi.

**Labels.** [Theorem] for 6.2 (gauge), 6.3 (frozen: only types I$_\pm$), 6.5 (conserved $\tau_d$) and 6.7 (the Sklyanin relation, with the involution as an algebraic corollary). List all four in the PR description as having had no referee pass.

**"Our presentation, not claimed new":**

- the block calculation and the pairing argument;
- the delay as a ratio of eigenvalues;
- the closed forms $\mathcal D_{\pm3}$ and the stored spin-3 charge;
- the hand proof that the window determines $r$;
- the refactorization;
- the $a_r$ component form.

## 7. Lectures 7–9: what the earlier lectures promise and what to read

**Promises to keep:**

- **Lecture 7.**
  - The exponentials of $\Phi$ and $\Psi$ stop commuting ($q$-Weyl) (L3 `sec-l3-lagrangian`).
  - The window keeps only half of the quantum group, the quantum counterpart of the chiral split (L3 `sec-l3-massless`).
  - Compare Konik–LeClair's quantum transmission matrix with the classical delay, and the odd window as a resonance at $\theta_d-\frac{i\pi}{2\lambda}$ (L4 `sec-l4-question`, `sec-l4-absorption`).
  - The lightest breather crosses with exactly the classical factor (L4 `sec-l4-waves`).
  - Wigner's phase versus the quantum transmission matrix (L4 `exr-l4-wigner`).
- **Lecture 8.** Wrap sine-Gordon windows around a cylinder and find Baxter's Q-operators (L1 `sec-l1-operators`). The chain of windows of 6.7 is a natural bridge.
- **Lecture 9.**
  - Repeat Bajnok–George's derivation of window crossing from boundary crossing-unitarity by folding (L1 `sec-l1-unitarity-crossing`).
  - Reconstruct the argument that fixes $\theta_b$ for the Lee–Yang window (L1 `sec-l1-lee-yang`).
  - Pressed against the boundary $\mathbb 1$, the window $\Phi(b)$ becomes the boundary $\Phi(b)$ (L1).
  - Why freezing $\Psi$ is not yet a derivation of the boundaries (L3 `exr-l3-boundary`).
  - The two ways out for other Toda theories: type-II windows, and fusing a window with its mirror image (L5 `sec-l5-others`). $a_r^{(1)}$ appears here.

**Book material.**

- Lecture 7: ch. 40 `sec-defect-qg` (Borel subalgebras, q-oscillators, transmission matrices as intertwiners, the two Borels and the two types) and ch. 41 (soliton transmission matrices, unitarity and crossing with the sign of $q$, the $a_2^{(2)}$ type-II window, bound states, semiclassics). Also ch. 37 `sec-sg-defect-spectrum` and the triangle relations of ch. 36.
- Lecture 8: ch. 40 `sec-universal-r-affine` and `sec-q-operators` (defect operators commute, the TQ relation, BLZ L-operators) and ch. 36 `sec-defect-circle`, `sec-topological-lines`, `sec-defect-finite-size`.
- Lecture 9: ch. 42 (the doubled theory, diagonal symmetric pairs, jump defects in the doubled picture, defect–boundary fusion, tensor K-matrices) and ch. 43 (breather and particle transmission, continuation, the defect spectrum). Also ch. 38 `sec-type-two-atft` and `sec-defect-folding`.
- The book's conventions for Part VIII are in `briefs/part8_defects_conventions.md`. Mind the traps in `AGENTS.md`, above all the sign of $q$ in the quantum lectures.

## 8. Follow-ups

1. PR #33 is merged, so Lectures 3 and 5 can now link the book's new ids `prp-dc-types`, `eq-dc-halves`, `eq-dc-lightcone` and `eq-dc-components`. In 03, sec-l3-types, the sentence that ends "(see its [defect Lagrangian](book:sec-defect-lagrangian))" should point to `prp-dc-types`. In 05, the sections `sec-l5-halves`, `sec-l5-types` and `sec-l5-cycle` can say the book now gives the same proof and the same relation between the two types. Small; offer it to the editor, or include it in the next lecture PR if the editor agrees.
2. One possible slip noticed in Delfino–Mussardo–Simonetti hep-th/9409076: eqs. (4.6) and (4.8) have $g/4m$ where their Lagrangian (4.1) gives $g/2m$. It was not raised in the book. Check before raising it.

## 9. Files of this handoff

- `code/lectures/defects_l6_lax.py`, `defects_l6_f1_holonomy.py`, `defects_l6_f2_charges.py`. All pass.
- `lecture-notes/defects/figures/l6_holonomy.svg` and `l6_charges.svg`.
- Rows in `code/lectures/README.md`, with the ids marked as planned.
