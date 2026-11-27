# Prior-art record: Chebyshev's bias for Ramanujan's tau function

Compiled 2026-11-24. Read-only literature survey; nothing here is a new claim.

Reading-status legend (applied honestly, per source):

- **[F]** full text read (statements and proofs, from the arXiv TeX source or the PDF).
- **[P]** partially read; the entry says exactly which parts.
- **[A]** abstract (or abstract page / search snippet) only.
- **[N]** not obtained.

Conventions. `lambda(p) = tau(p) p^{-11/2} = 2 cos(theta_p)`, so `tau(p)/p^6 = lambda(p)/sqrt(p)`. `L(s, Delta) = sum tau(n) n^{-s}`; the normalised function is `L(s + 11/2, Delta)` with centre `s = 1/2`. Quotes are verbatim from the TeX source where one was available (LaTeX markup kept); `[...]` marks an omission of mine. Where a statement in this file is my own deduction or remark and not in a source, it is labelled **(compiler's remark)**.

Downloads (TeX sources, PDFs) were kept outside the repository, in the session scratch directory; none is committed.

---

## 0. Access log (what could and could not be fetched)

| Source | Result |
|---|---|
| arXiv abs pages, `arxiv.org/e-print/` TeX sources, PDFs for 2203.12791, 2509.04187 (v1 and v2), 2405.01512, 2512.23302 (v3), 2607.28931 (v2), 2607.23150, 2411.07436 (v3), 1706.06394 (v4), 2209.12832, 2203.12266, 2206.02612, math/9411213, 2002.04717 | all fetched |
| Published version of Sheth (Cambridge Core, open access, doi 10.1017/S0305004125000234) | fetched (PDF) |
| Sarnak's letter, `publications.ias.edu/sites/default/files/MazurLtrMay08.PDF` | HTTP 403 (Cloudflare challenge page). Obtained instead from `web.math.princeton.edu/sarnak/MazurLtrMay08.PDF` (6-page PDF, 11 numbered pages of letter incl. the January 2008 addendum) |
| Rubinstein-Sarnak 1994, Project Euclid | blocked (Incapsula bot page). Obtained from `www.math.uwaterloo.ca/~mrubinst/publications/Chebyshev.pdf` (25 pages; Type-3 fonts, so text extraction fails; pages were read as images) |
| web.archive.org | connection reset |
| LMFDB `www.lmfdb.org` HTML pages via curl | Google reCAPTCHA challenge page on every L-function URL tried |
| LMFDB JSON API | one call succeeded (first 100 degree-3, conductor-1 records: all GL(3) Maass forms); subsequent calls returned the reCAPTCHA page. I did not try to get around the challenge |
| LMFDB via the web-fetch tool (which returns a model-written summary of the page, not raw HTML) | worked for four pages, then returned "Checking your browser before accessing www.lmfdb.org"; details in section 9.4 |
| Published version of Koyama-Kurokawa (Proc. Japan Acad.) | not fetched; arXiv v2 read instead |
| Yoshida 1988 (J. Ramanujan Math. Soc. 3) | not obtained |

---

## 1. Koyama-Kurokawa, tau and DRH  **[F]** (arXiv v2)

S. Koyama, N. Kurokawa, "Chebyshev's Bias for Ramanujan's tau-function via the Deep Riemann Hypothesis", arXiv:2203.12791 (v1 24 Mar 2022, v2 25 Mar 2022). Journal data as given in the bibliographies of Sheth and Hayani: Proc. Japan Acad. Ser. A Math. Sci. 98 (2022), no. 6, 35-39, doi 10.3792/pjaa.98.007, MR4432981. The arXiv abs page lists no journal reference. I read the arXiv v2 TeX source in full (it is 7 KB); the published version was not compared.

Hypothesis used: the Deep Riemann Hypothesis (convergence of the Euler product at the centre), which the paper itself notes is stronger than RH. Verbatim:

> \textbf{Deep Riemann Hypothesis (DRH):} The Euler product \eqref{EP} converges at $s=\f12$ with
> $$\lim_{x \to \infty}\prod_{p \le x} \det\l(1-M(p)p^{-\f12}\r)^{-1} =\sqrt{2}^{\delta(M)} L\l(\tfrac12, M\r),$$
> where $\delta(M)=\mathrm{ord}_{s=1}L(s,M^2)$ with ord$_{s=1}$ signifying the order of the pole at $s=1$.

Standing assumptions (Section 1): `L(s,M)` entire, functional equation `s <-> 1-s`, and "Moreover, $L(\f12,M)\ne0$ is assumed for simplicity." For Delta: "$\delta(M)=-1$" (Example 2).

> **Theorem 1.** Suppose DRH for $L(s,M)$, then
> $$\lim_{x\to\infty}\f{\displaystyle\sum\limits_{p\le x}\f{\tr(M(p))}{\sqrt p}}{\log\log x}=-\f{\delta(M)}2.$$

> **Theorem 2.** Assume the DRH for $L(s+\f{11}2,\Delta)$. The following holds.
> $$\sum_{\genfrac{}{}{0pt}{1}{p\le x}{\text{prime}}}\f{\tau(p)}{p^6}\sim\f12\log\log x\quad(x\to\infty).$$

> **Corollary 1.** Assume the DRH for $L(s+\f{11}2,\Delta)$. The natural density of the set
> $$A=\l\{x>0\ \l|\ \sum_{p\le x}\f{\tau(p)}{p^6}>0\r\}\r.$$
> exists and is equal to 1. More precisely,
> $$\lim_{X\to\infty}\f1X\int_{t\in A\cap[2,X]}dt=1.$$

Proof of Theorem 1 is the three-line decomposition `I + II + III = log(partial Euler product) = O(1)` with `II(x) ~ (delta(M)/2) log log x` (a generalised Mertens theorem cited from Kaneko-Koyama-Kurokawa) and `III = O(1)`.

What the paper says about Sarnak (verbatim): "Sarnak \cite{S} also reached a similar prediction under the assumption of the Generalized Riemann Hypothesis as well as the Grand Simplicity Hypothesis for $L(s,\Delta)$ [...] He has pointed out that the sum $S(x)=\sum_{p\le x}\f{\tau(p)}{p^\f{11}2}$ has a bias to being positive, in the sense that the mean of the measure $\mu$ [...] is equal to 1."

Not in the paper: any statement under GRH alone; any explicit `x_0`; any error term; any converse.

Restatement elsewhere **[P]**: Aoki-Koyama, "Chebyshev's bias against splitting and principal primes in global fields", arXiv:2203.12266 (J. Number Theory 245 (2023) 233-262), Section "Perspectives", restates this as "Theorem (Koyama-Kurokawa)": "Assume DRH(A) for $L(s+\f{11}2,\Delta)$. Then the sequence $a(p)=\tau(p)p^{-11/2}$ has a Chebyshev bias to being positive. More precisely, it holds that $\sum_{p\le x}\f{\tau(p)}{p^6}=\f12\log\log x+c+o(1)\quad(x\to\infty)$ with some constant $c$." I read only that subsection of Aoki-Koyama.

---

## 2. Koyama-Sheth, modular forms  **[F]** for Sections 1-2 and all theorem statements; **[P]** for Section 3 proofs and the appendix (skimmed)

S. Koyama, A. Sheth, "Chebyshev's bias for modular forms", arXiv:2509.04187, v1 4 Sep 2025, v2 15 Jan 2026 ("Updated version; minor typos corrected"). No journal reference on the abs page. I diffed the v1 and v2 TeX sources: the differences are an address change, `\leq` to `\leqslant`, a label rename, one added sentence of references, and typo fixes. No statement changed.

Setting: `f` a non-CM cuspidal newform on `Gamma_0(N)`, trivial nebentypus, normalised coefficients `a_f(p) = 2 cos theta_p`.

> **Theorem 1.5.** Let $f \in S_k(\Gamma_0(N))$ be a cuspidal newform. Assume Conjecture \ref{DRH} for $L(s, f)$ and let $m(f)=\ord_{s=1/2} L(s, f)$. Then there exists a constant $c_f$ such that
> $$\sum_{p \leqslant x} \frac{a_f(p) }{\sqrt p} = \left( \frac{1}{2}- m(f) \right) \log \log x +c_f +o(1).$$
> In particular, in the sense of Definition \ref{Chebyshev}, we conclude that (1) The sequence $(a_f(p))_p$ has a Chebyshev's bias towards being positive if $m(f)=0$. (2) The sequence $(a_f(p))_p$ has a Chebyshev's bias towards being negative if $m(f) \geq 1$.

"Conjecture \ref{DRH}" is the Kaneko-Koyama-Kurokawa Deep Riemann Hypothesis (their Conjecture 1.1). So the main theorem is under DRH, not GRH.

The one sentence about GRH alone (in the paragraph after Theorem 1.5), verbatim:

> we remark that using \cite[Theorem 4.1]{Sheth24B}, the asymptotic in Theorem \ref{mainthm6} holds as $x \to \infty$ outside a set of finite logarithmic measure assuming only GRH, and so we can in fact obtain a statement which captures the bias assuming only GRH.

And in the introduction: "the second author showed in \cite{Sheth24B} that the GRH implies Conjecture \ref{DRH} as $x \to \infty$ outside a set of finite logarithmic measure."

Tau appears as Example 1.4 (Koyama-Kurokawa's theorem, under DRH) and in:

> **Theorem 3.9.** The following hold true. (1) For $k \in \{12, 16, 20, 18, 22, 26\}$, $\Lambda(\Delta_k, s)$ satisfies the super-positivity property. (2) If $k=12, 16$ or $20$, then $m(\Delta_k)=0$. (3) If $k=18, 22$ or $26$, then $m(\Delta_k)=1$.

> **Corollary 3.10.** Assume Conjecture \ref{DRH} for $L(s, \Delta_k)$. We have that
> $$\sum_{p \leqslant x} \frac{\tau_{k}(p)}{p^{\frac{k}{2}}} = \begin{cases} \dfrac{1}{2} \log \log x+c_k+o(1) & \text{if $k = 12, 16, 20 $} \\ -\dfrac{1}{2} \log \log x+c_k+o(1) &\text{if $k = 18, 22, 26$} \end{cases}$$
> for some constant $c_k$. In particular [...] the sequence $(\tau_k(p) p^{-\frac{k-1}{2}})_p$ has Chebyshev's bias towards being positive if $k \in \{12, 16, 20\}$ and a Chebyshev's bias towards being negative if $k \in \{18, 22, 26\}$.

Here `Delta_k = Delta * E_{k-12}`. The proof that `m(Delta) = 0` is unconditional: `Lambda(Delta, 1/2)` is an integral of the positive function `Delta(iy)`.

Proposition 2.2: `R(pi_f) = 1` for every cuspidal newform (second-moment L-function has a simple zero at `s = 1`: `Sym^2` is holomorphic and non-vanishing at 1 by Shahidi, `wedge^2 = zeta`).

Not in the paper: no all-`x` statement under GRH; no explicit constants or thresholds ("explicit" in this paper means "explicit value of `m(f)`"); no converse; nothing about smoothed sums.

---

## 3. Sheth, Euler products at the centre  **[F]** (arXiv v1 TeX in full; published PDF: theorem statements and Section 4 compared)

A. Sheth, "Euler Products at the Centre and Applications to Chebyshev's Bias", arXiv:2405.01512 (single version, 2 May 2024). Published: Math. Proc. Camb. Phil. Soc. 179 (2025), 331-349, doi 10.1017/S0305004125000234, received 8 May 2024, revised and accepted 30 December 2024, first published online 27 May 2025, open access.

**Numbering differs between versions.** Koyama-Sheth and Hayani cite the published numbering.

| Content | arXiv v1 | published |
|---|---|---|
| GRH implies DRH off a finite-log-measure set (Theorem A) | Theorem 3.8 | Theorem 3.9 |
| `psi(x,pi) << sqrt(x) (log log x)^2` off a finite-log-measure set | Theorem 3.4 | Theorem 3.5 |
| generalised Mertens lemma | Lemma 3.5 | Lemma 3.6 |
| Chebyshev's bias for Satake parameters (Theorem B) | Theorem 4.2 | **Theorem 4.1** |
| mod 4 race (Corollary A) | Corollary 4.3 | Corollary 4.2 |
| **tau** (Corollary B) | Corollary 4.4 | **Corollary 4.3** |
| elliptic curves (Corollary C) | Corollary 4.5 | Corollary 4.4 |

Which L-functions: irreducible cuspidal automorphic representations `pi` of `GL_n(A_Q)` with `L(s,pi)` entire. Hypotheses: the Riemann Hypothesis for `L(s,pi)` and the Ramanujan-Petersson conjecture for `pi`. No linear independence, no DRH. Exceptional set: finite logarithmic measure, with

> **Definition 1.4.** Let $S \subseteq \mathbb R_{\geq 2}$ be a measurable subset of the real numbers. The logarithmic measure of $S$ is defined to be $\mu^{\times}(S)= \int_{S} \frac{dt}{t}.$

Verbatim statements (arXiv TeX):

> **Theorem A** (arXiv Theorem 3.8). Let $\pi$ be an irreducible cuspidal automorphic representation of $\GL_n(\mathbb A_\mathbb Q)$ such that $L(s, \pi)$ is entire and let $m = \ord \limits _{s = \frac{1}{2}} L(s, \pi)$. Assume the Ramanujan--Petersson Conjecture and the Riemann Hypothesis for $L(s, \pi)$. Then there exists a subset $S \subseteq \mathbb R_{\geq 2}$ of finite logarithmic measure such that for all $x \not \in S$,
> $$(\log x)^m \cdot  \prod_{p \leq x}  \prod_{j=1}^n (1-\alpha_{j, p}p^{-\frac{1}{2} })^{-1}   \sim  \frac{\sqrt{2}^{ \nu (\pi) }}{e^{m \gamma} m!} \cdot L^{(m)} \left(\frac{1}{2}, \pi \right).$$

> **Theorem B** (arXiv Theorem 4.2, "Chebyshev's bias for Satake parameters"). Let $\pi$ be an irreducible cuspidal automorphic representation of $\GL_n(\mathbb A_\mathbb Q)$ such that $L(s, \pi)$ is entire and let $m= \ord_{s=1/2} L(s, \pi)$. Assume the Ramanujan--Petersson Conjecture and the Riemann Hypothesis for $L(s, \pi)$. Then there exists a constant $c_\pi$ such that
> $$\sum_{p \leq x} \frac{\alpha_{1, p}+ \cdots +\alpha_{n, p}}{\sqrt p}= \left( \frac{R(\pi)}{2}-m \right) \log \log x+ c_\pi+ o(1)$$
> for all $x$ outside a set of finite logarithmic measure, where $R(\pi) =\ord \limits _{s=1} L_2(s, \pi)$.

In the published version the left-hand side of Theorem B / Theorem 4.1 carries a real part, `Re( sum ... )`, and the hypotheses read "the Generalised Riemann Hypothesis and the Ramanujan-Petersson Conjecture for L(s, pi)". (The published PDF's text layer drops some symbols, so I quote math from the arXiv TeX and record only the wording differences.)

> **Corollary B** (arXiv Corollary 4.4; published Corollary 4.3). Let $\tau(n)$ denote Ramanujan's tau function. Assume the Riemann Hypothesis for $L(s, \Delta)$. Then there exists a constant $c$ such that
> $$\sum_{p \leq x} \frac{\tau(p)}{p^6} = \frac{1}{2} \log \log x + c+o(1)$$
> for all $x$ outside a set of finite logarithmic measure. In particular, in the sense of Definition \ref{Chebyshev}, the sequence $\tau(p) p^{-\frac{11}{2}}$ has a Chebyshev bias towards being positive.

Proof of the corollary (verbatim, arXiv): "By the work of Deligne \cite{Deligne1974} we have that $|\tau(p)| \leq 2 p^{11/2}$ [...] It is known that $L \left(s, \pi_\Delta \right) \neq 0$ and $R(\pi_\Delta)=1$ (see \cite{KK22}). Thus, the asymptotic follows from Theorem \ref{Satake}". So for Delta the Ramanujan-Petersson hypothesis is Deligne's theorem and the only unproved hypothesis is RH for `L(s, Delta)`.

Where the exceptional set comes from:

> **Theorem 3.4** (arXiv; published 3.5). Let $\pi$ be an irreducible cuspidal representation of $\GL_n(\mathbb A_\mathbb Q)$ such that $L(s, \pi)$ is entire. Assume the Riemann Hypothesis for $L(s, \pi)$. Then there exists a set $S$ of finite logarithmic measure such that for all $x \not \in S$, $\psi(x, \pi) \ll \sqrt x (\log \log x)^2.$

Its proof in the paper is a citation: "The analog of this result for the prime counting function $\psi(x)$ was first proven by Gallagher \cite{Gallagher1980} and was later generalized in the setting of the theorem by Qu (\cite[Theorem 1.1]{Qu2007}); similar arguments have also been given in \cite[Theorem 2]{Koyama2016} and \cite[Theorem 3.4]{Sheth2023}." The explicit-formula asymptotic (arXiv Theorem 3.2) is likewise proved by reference: "These integrals can now be analysed in an analogous way as in \cite[Theorem 2.3]{Sheth2023}". I did not read Gallagher 1980, Qu 2007 or Sheth 2023 (arXiv:2312.05236).

The author's own novelty claim (verbatim): "The key difference is that the statements in \textit{op.cit.} assume Conjecture \ref{DRH}, but our statements only assume the Generalized Riemann Hypothesis, with the caveat that our asymptotics hold for all $x$ outside an exceptional set. Indeed, as far as we are aware, this is the first instance in the literature where an explanation of the above mentioned phenomenona can be given only assuming the Generalized Riemann Hypothesis".

**Answer to the question asked.** Yes: under RH for `L(s, Delta)` alone, Corollary B gives `sum_{p<=x} tau(p)/p^6 = (1/2) log log x + c + o(1)` for all `x` outside a set `S` with `int_S dt/t < infinity`. Consequences:

- **(compiler's remark, immediate)** the set `{x : sum_{p<=x} tau(p)/p^6 <= 0}` is contained in `S` union a bounded set, hence has finite logarithmic measure, hence logarithmic density 0.
- Natural density: the paper does not mention natural density. The step "finite logarithmic measure implies natural density 0" is written out in Hayani (section 4 below, his display (1.5)). So "natural density 1 under GRH alone" for tau is a two-line corollary of Sheth plus that remark, but I did not find that sentence written for tau anywhere.

What the paper does not contain: no rate for the `o(1)`; no value or bound for `c`; no description of `S` beyond finiteness of its logarithmic measure; no all-`x` statement; no statement for smoothed sums.

**A converse that is in the TeX source but not in the paper.** The arXiv TeX contains, inside `\begin{comment} ... \end{comment}` blocks (so it is absent from the arXiv PDF and from the published version), a draft theorem: "Assume the Ramanujan conjecture for $L(s, \pi)$. Suppose that $(\log x)^m \prod_{p \leq x}  \prod_{j=1}^n (1-\alpha_{j, p}p^{-\frac{1}{2} })^{-1}   \asymp 1.$ Then the Riemann Hypothesis for $L(s, \pi)$ holds." with a proof sketch "Using the method of Goldfeld". This is commented-out draft material, not a published result, and should not be cited as one. (Conrad 2005 and Goldfeld 1982, which I did not read, are the cited sources for "convergence of the Euler product implies GRH".)

---

## 4. Hayani, Chebyshev's bias without linear independence  **[F]** (v3 TeX, including the appendix)

M. Hayani, "Chebyshev's bias without linear independence", arXiv:2512.23302, v1 29 Dec 2025, v2 5 Jan 2026, v3 25 Jun 2026. No journal reference on the abs page. I read v3 only.

Scope: Dirichlet characters modulo a fixed `q`. Notation: `pi_{1/2}(x;q,a) = sum_{p<=x, p=a mod q} p^{-1/2}`, `M(q;a,b) = (r(a)-r(b)+2 sum_{chi != chi_0}(chi(a)-chi(b)) m_chi) / (2 phi(q))`, `r(a)` the number of square roots of `a` mod `q`, `m_chi = ord_{s=1/2} L(s,chi)`.

> **Theorem 1.1.** Let $q\ge 3$ and assume \text{GRH} for all non-principal Dirichlet characters modulo $q$. Let $\varepsilon>0$ and let $a,b$ be distinct invertible residue classes modulo $q$. Then there exists a (unique) constant $C$ depending on $q,\, a,$ and $b$, such that the natural density of the set
> \[\left\{x\ge 2\, :\, \left| \pi_{1/2}(x;q,a)-\pi_{1/2}(x;q,b) +\M(q;a,b)\log \log x-C\right|\le \frac{(\log \log x)^{3+\varepsilon}}{\log x}\right\}\]
> exists and equals $1$. Furthermore, there exists $K=K(a,b)>1$ such that the logarithmic density of the set
> \[\left\{x\ge 2\, :\, \left| \pi_{1/2}(x;q,a)-\pi_{1/2}(x;q,b) +\M(q;a,b)\log \log x-C\right|\le \frac{K \log \log x}{\log x}\right\}\]
> exists and equals $1$.

> **Corollary 1.2.** Assume GRH for the non-principal character modulo $4$. Then, the natural density of the set $ \mathcal{P}_{1/2}(4;3,1)=\{ x\ge 2\, :\, \pi_{1/2}(x;4,3)>\pi_{1/2}(x;4,1)\}$ exists and equals $1$.

> **Theorem 1.3.** Assume GRH for a non-principal Dirichlet character $\chi$ modulo $q$. Then, there exists $\ell_\chi\ne 0$ such that, for all $\varepsilon>0$, the natural density of the set
> $$\mathcal{E}_\chi:=\left\{x\ge 2\, :\, \left|  (\log x)^{m_\chi} \prod_{p\le x}\left(1-\frac{\chi(p)}{p^{1/2}}\right)^{-1}-\ell_\chi  \right|\le \frac{(\log \log x)^{3+\varepsilon}}{\log x} \right\}$$
> exists and equals $1$.

> **Theorem 1.4.** Let $q\ge 3$ and assume \text{GRH} for all non-principal Dirichlet characters modulo $q$. Let $a,b$ be distinct invertible residue classes modulo $q$. Then
> \[\frac{1}{x}\int_{ 2}^x \bigl(\pi_{1/2}(u;q,a)-\pi_{1/2}(u;q,b)\bigr)\mathrm{d}u=-\M(q;a,b)\log \log x+C+O\left(\frac{\log \log x}{\log x}\right)\, .\]

Theorem 1.4 is an **every-`x`** statement under GRH alone for a Cesaro mean. **(compiler's remark)** Since `(1/x) int_2^x sum_{p<=u} c_p du = sum_{p<=x} c_p (1 - p/x)`, the left side of Theorem 1.4 is exactly a `(1 - p/x)`-smoothed sum of `(1_{p=a} - 1_{p=b}) p^{-1/2}`. This is the Dirichlet-character analogue of item 10(a) in asymptotic form: eventual one-signedness under GRH when `M(q;a,b) != 0`, with a non-explicit `O`-constant and no explicit `x_0`.

Does it cover modular forms or tau? **No.** Tau and modular forms are not treated. The only extension the paper mentions is: "Although our arguments extend straightforwardly to the context of global fields, we restrict our focus to the classical case of primes in arithmetic progressions". The proof uses the explicit formula for `pi(x;t)` via Rubinstein-Sarnak Lemma 2.1 and a `2k`-th moment bound (Lemma 2.1, `(1/Y) int |Delta(y;t)|^{2k} dy <= (C k)^{4k}`, proved in the appendix after Puchta) that is written for zeros of Dirichlet L-functions.

Relation to Sheth, from the abstract (verbatim): "The same type of results presented here were independently proved by Arshay Sheth (2025) in the general context of automorphic forms, which implies our main asymptotic. While the spirit of the proofs is similar, Sheth develops an explicit formula for the partial Euler product and uses a result due to Gallagher (1980) to prove that some estimates hold outside a set of finite logarithmic measure. We share this independent work because it provides a completely self-contained and elementary proof relying only on the usual explicit formula, and it yields explicit error terms rather than an implicit $o(1)$ asymptotic."

"Explicit error terms" means the shape `(log log x)^{3+eps}/log x`; the constants `C`, `K`, `Y_0` and the implied constants are not numerically explicit.

The measure-theoretic remark needed in section 3 above (verbatim, his (1.5)): "Indeed, any Borel set $A\subset (0,\infty)$ satisfying $\int_0^\infty \mathds{1}_A(e^y)\mathrm{d}y<\infty$ has zero natural density. To see this, note that for any $\mathcal{L}>1$,
$$0\le \limsup_{x\to \infty}\frac{1}{x} \int_1^x \mathds{1}_A(u)\mathrm{d}u\le \limsup_{x\to \infty}\int_{\mathcal{L}}^x \mathds{1}_A(u) \frac{\mathrm{d}u}{u}=\int_{\log \mathcal{L}}^\infty \mathds{1}_{A}\bigl(e^y\bigr)\mathrm{d}y\, .$$
The result follows by letting $\mathcal{L}\to \infty$."

His description of Suzuki (useful for section 6): "The author proved that GRH holds for the non-principal Dirichlet character modulo $4$ if and only if [the `log p log(x/p)`-weighted difference] $\sim-(\log x)^2 \frac{\sqrt{x}}{4}$".

---

## 5. The two July 2026 papers

### 5.1 Koyama, arXiv:2607.28931  **[P]** (abstract, introduction and all statements of Section 1 read; proofs not read)

S. Koyama, "The Fine-Structure Hierarchy of Prime Biases and the Universal Dominance of $-1 \pmod N$", arXiv:2607.28931, v1 31 Jul 2026, v2 29 Aug 2026. Comments: "10 pages, 4 tables. This is a revised and expanded version of the author's previous manuscript 'A Hidden Hierarchy of Chebyshev's Bias and the Dominance of -1 (mod N)' (submitted April 26, 2026)".

Scope: Dirichlet L-functions mod `N`; Gaussian-mollified Weil explicit formula; hypothesis DRH.

> **Theorem 1.4.** Assume DRH for all Dirichlet $L$-functions $L(s,\chi)$ modulo $N$. For any fixed $T>0$, as $x \to \infty$, the spectrally normalized, mollified fine-structure prime bias difference $\widetilde{S}_T(x, \chi_{1,a})$ satisfies:
> \[\widetilde{S}_T(x, \chi_{1,a}) = C_N \cdot \log L(1, \chi_{1,a}) + \mathcal{O}\left(\frac{\log x}{\sqrt x}\right),\]
> where $C_N > 0$ is an explicit positive constant depending solely on the modulus $N$ [...]

> **Corollary 1.6** ("Universal Determination of Multi-Way Ranking Hierarchy"), item 1: "the inequality $S_T(x,a) > S_T(x,1)$ holds for all sufficiently large $x$ if and only if $\log L(1, \chi_{1,a}) < 0$."

(Theorem and corollary numbers are those of the v2 PDF. In the TeX source of v2 the displayed definition of the normalised sum `\widetilde S_T(x,a)` is empty apart from a leftover author note, so the normalisation is not stated there.)

Tau / modular forms / automorphic: a case-insensitive search of the TeX source for `modular`, `ramanujan`, `\tau`, `cusp form`, `automorphic` returns **0 matches**. GRH-only: no; the paper states that GRH is insufficient for its argument ("Remark [Sharpness of the GRH Bound and Necessity of DRH]"). All-`x`: only "for all sufficiently large `x`", under DRH, for a mollified statistic. Explicit thresholds: none.

### 5.2 Hayani, arXiv:2607.23150  **[P]** (abstract, Section 1 statements, the weighted-races subsection; group theory not read)

M. Hayani, "On the role of higher roots in prime ideal races", arXiv:2607.23150, 25 Jul 2026, report number MPIM-Bonn-2026.

Scope: prime ideal races in Galois extensions `L/K` of number fields (Artin L-functions). Relevant statements:

> "More recently, Sheth~\cite{Sheth} weakened this assumption, proving that under $\mathrm{GRH}(L)$, the equality \eqref{DRHequality} holds outside an exceptional set of finite logarithmic measure."

> **Theorem 1.6.** Let $L/K$ be a Galois extension of number fields with Galois group $G$ and let $t\, \colon\,  G\to \R$ be a non-zero class function such that $\iota(t)<\infty$ and $\gamma:=\gamma(t)<2\iota(t)$. If $\mathrm{GRH}(L)$ holds, then
> $$\pi_\alpha(x;L/K;t)=\begin{cases} K_{\alpha,t}\frac{x^{(1-\alpha \gamma)/\gamma}}{\log x}+O\left(x^{(1-\alpha \gamma)/\gamma}(\log x)^{-2}\right) &\text{if}\ \  \alpha<1/\gamma\\ K'_t \log \log x+O(1) &\text{if}\ \alpha=1/\gamma\\ O(1)  &\text{if}\ \ \alpha > 1/\gamma \end{cases}$$
> where $K_{\alpha,t}$ and  $K'_t$ are non-zero constants having the same sign as $\mu(\gamma)\langle t,r_\gamma\rangle_G$. If, moreover, $\gamma(t)=\iota(t)$, then~\eqref{pialpha-asympt} holds unconditionally.

> **Theorem 1.7.** Let $p\ge2$ be a prime. There exist infinitely many Galois extensions of number fields $L/K$ with Galois group $G$ and conjugacy classes $C_1,C_2\in G^\#$, such that for all $x\ge2$ we have [the same three-case asymptotic for $\frac{1}{|C_1|}\pi_{\alpha}(x;L/K;C_1)-\frac{1}{|C_2|}\pi_{\alpha}(x;L/K;C_2)$ with exponent $1/p$] where $K_{\alpha,C_1,C_2}, K'_{C_1,C_2}$ are positive constants.

These are all-`x` asymptotics (some unconditional) for special Artin-type races where a class function is orthogonal to low-order root-counting functions. They are of the "extreme bias" type. Nothing about tau, modular forms or cuspidal `GL_2`: a search of the source for `modular`, `ramanujan`, `cusp`, `automorphic` finds no relevant hit. No numerically explicit thresholds.

### 5.3 Other 2026 items found by the searches  **[A]**

- M. Aoki, "Distribution of Primes in Galois Extensions of Global Fields with Explicit Constant Terms", arXiv:2611.20592 (16 Nov 2026). Abstract: a refined asymptotic for weighted prime counting functions "is equivalent to the full version of DRH, namely the convergence of the Euler products on the critical line (DRH(A)) and the determination of their limiting values (DRH(B)). Here, the constant term is given explicitly". Artin L-functions; DRH; "explicit" refers to the constant term `c`. Not tau.
- S. Han, "On variants of Polya's conjecture", arXiv:2608.27130 (27 Aug 2026). Abstract: for `f(x) = - sum_{n<=x} lambda(n) log n n^{-1/2} log(x/n)` (Liouville), "if $f(x)$ is non-negative for all sufficiently large $x$, then the Riemann Hypothesis holds"; the converse needs RH plus simple zeros plus a convergence hypothesis. Zeta only.
- N. K. Sinha, "Chebyshev Bias for Largest Prime Factors", arXiv:2608.06435 (6 Aug 2026). Abstract: under RH for `zeta` and `L(s, chi_{-4})`, "$D_k(x)<0$ for all sufficiently large $x$". Not tau.
- A. Sheth, W. Voo, "Super-positivity via Dirichlet theta functions", arXiv:2611.12152 (9 Nov 2026). Quadratic characters. Not tau.

---

## 6. Suzuki, variants of Chebyshev's conjecture

### 6.1 arXiv:2411.07436  **[P]** (v3 TeX: Section 1 and every theorem statement, including Theorems 8 and 9 in Section 4, read in full; proofs not read)

M. Suzuki, "On variants of Chebyshev's conjecture", arXiv:2411.07436, v1 11 Nov 2024, v3 10 Nov 2025. The Ramanujan Journal 68, article number 95 (2025), doi 10.1007/s11139-025-01238-9.

Scope: the Riemann zeta function and Dirichlet L-functions only. A search of the source for `automorphic`, `modular`, `ramanujan`, `cusp`, `elliptic`, "Selberg class" finds only the Gauss sum `\tau(\chi)` and one bibliography title. **No statement about tau or automorphic L-functions.**

Statements of the form "one sign for all `x >= x_0` iff (G)RH":

> **Theorem 1.** The following three statements are equivalent: (1) The RH holds. (2) There exists an $x_0 \geq 2$ such that
> $$\int_{0}^{x} (\psi_{1/2}(y) - 2\sqrt{y}\,) \frac{dy}{y} = \sum_{n \leq x} \frac{\Lambda(n)}{\sqrt{n}}\log\frac{x}{n}  - 4 \sqrt{x} ~\leq~  0$$
> holds for all $x \geq x_0$. (3) It holds that
> $$\lim_{x \to \infty}\left[\sum_{n \leq x} \frac{\Lambda(n)}{\sqrt{n}}\left(1-\frac{\log n}{\log x} \right)- \frac{4\sqrt{x}}{\log x}\right] = -\frac{\zeta'}{\zeta}\left(\frac{1}{2}\right).$$

> **Theorem 3.** The following three statements are equivalent: (1) The GRH for $L(s,\chi_4)$ holds. (2) There exists an $x_0 \geq 2$ such that
> $$\int_{0}^{x} (\psi_{1/2}(y;4,1)-\psi_{1/2}(y;4,3)) \frac{dy}{y} = \sum_{n \leq x} \frac{\Lambda(n)\chi_4(n)}{\sqrt{n}}\log\frac{x}{n} ~\leq~  0$$
> holds for all $x \geq x_0$. (3) It holds that
> $$\lim_{x \to \infty}\sum_{n \leq x} \frac{\Lambda(n)\chi_4(n)}{\sqrt{n}}\left(1- \frac{\log n}{\log x}\right) = -\frac{L'}{L}\left(\frac{1}{2},\chi_4\right).$$

> **Theorem 4.** The GRH for $L(s,\chi_{4})$ holds if and only if
> $$\lim_{x \to \infty} \sum_{2<p \leq x} (-1)^{(p-1)/2} \log p \cdot\sqrt{\frac{x}{p}}\,\log\frac{x}{p} = -\infty.$$

Theorem 4 is a sum over **primes only**. The general versions are Theorem 8 (any non-principal `chi`; sums over prime powers) and Theorem 9 (real `chi`; sums over primes). From Theorem 8, verbatim:

> Suppose that there exists $x_0 \geq 2$ such that $\sum_{n \leq x} \frac{\Lambda(n)}{\sqrt{n}}\Re(\chi(n))\log\frac{x}{n}$ has a constant sign for all $x \geq x_0$, and $L(\sigma)\not=0$ for $\sigma>\beta$ for some $1/2 \leq \beta <1$. Then, $L(s)\not=0$ in the right-half plane $\Re(s)>\beta$. [...] Suppose that $L(1/2)\not=0$. Then, we have
> $$\lim_{x \to \infty}\sum_{n \leq x} \frac{\Lambda(n)\chi(n)}{\sqrt{n}}\left(1- \frac{\log n}{\log x}\right) = -\frac{L'}{L}\left(\frac{1}{2}\right)$$
> if and only if the GRH for $L(s)$ holds. If $\chi$ is primitive, then equation \eqref{EQ_402} implies that the sum \eqref{EQ_401} has a constant sign for all sufficiently large $x>0$.

From Theorem 9 (real `chi`, `L(sigma,chi) != 0` for `sigma > beta >= 1/2`), verbatim: "Suppose that $\sum_{p \leq x} \chi(p) \log p \cdot \sqrt{\frac{x}{p}}\,\log\frac{x}{p}  \leq 0$ for all sufficiently large $x>1$. Then, $L(s,\chi) \not=0$ in the right-half plane $\Re(s)>\beta$. In particular, if $\beta=1/2$, the GRH for $L(s,\chi)$ holds. Conversely, suppose that the GRH for $L(s,\chi)$ holds and $L(1/2,\chi) \ne 0$. Then, the asymptotic formula $\sum_{p \leq x} \chi(p) \log p \cdot \sqrt{\frac{x}{p}}\,\log\frac{x}{p}  = -\frac{\sqrt{x}}{4}(\log x)^2\left(1+O\left(\frac{1}{\log x}\right) \right)$ holds as $x \to \infty$."

Other theorems: Theorem 2 and Theorem 5 (zeta, sums up to `x e^2`: sign constancy implies RH; the converse needs RH plus `sum_rho x^rho/rho = o(sqrt(x) log x)`); Theorem 6 (residue class `1 mod q`); Theorem 7 (unconditional negativity on average over `3 <= q <= Q`, `Q >= x`). Corollary 1: RH iff the non-archimedean part `g_0` of the zeta screw function is eventually of one sign.

All thresholds `x_0` are existential; none is explicit. The paper notes "Considering \cite{LS94}, we are interested in the logarithmic density of $x$ for which $\psi_{1/2}(x)-2\sqrt{x} \leq 0$ holds, but we do not study that topic in this paper."

Prior results as described in Suzuki's introduction (I did not read the originals): Johnston (arXiv:2201.06184, **[A]**: "the Riemann hypothesis is equivalent to the condition $\int_{2}^x\left(\pi(t)-\text{li}(t)\right)\mathrm{d}t<0$ for all $x>2$" - an all-`x`, explicit-threshold equivalence, for zeta); Hardy-Littlewood and Landau (Chebyshev's `e^{-p/x}`-weighted conjecture is equivalent to GRH for `L(s, chi_4)`); Knapowski-Turan; Fujii; Platt-Trudgian. All concern zeta or `chi_4`.

### 6.2 Suzuki, screw functions for the extended Selberg class, arXiv:2209.12832  **[P]** (abstract, Section 1, the Selberg-class axioms, statement of Theorem 4.1; proofs and Section 5 not read)

M. Suzuki, "Screw functions of Dirichlet series in the extended Selberg class", arXiv:2209.12832 (v1 14 Sep 2022). International Journal of Number Theory 21 (2025), no. 8, 1815-1823, doi 10.1142/S1793042125500885 (journal data from the arXiv abs page).

This is the closest thing I found to a "GRH iff eventual sign" statement that covers `L(s, Delta)`. The introduction says: "The Riemann zeta-function, Dirichlet $L$-functions, and $L$-functions of holomorphic cusp forms are typical members of $\mathcal{S}$". Tau and Delta are not singled out.

> **Theorem 1.1.** Let $F$ be a member of $\mathcal{S}^{\sharp}$ and let $g_F(t)$ be the screw function defined by \eqref{eq_101}. We assume that $F(s)$ has no real zeros except for the possible zero at $s=1/2$. Then, the GRH for $F$ is true if and only if $\Re(-g_F(t)) \geq 0$ for all $t \geq t_0$ for some $t_0 \geq 0$.

with `g_F(t) := - i B_F t - (m_0/2) t^2 + sum_{gamma != 0} m_gamma (e^{-i gamma t} - 1)/gamma^2` (sum over zeros of `F(1/2 - iz)`), and the zero-free formula:

> **Theorem 4.1.** Let $F \in \mathcal{S}^{\sharp\flat}$ [...] define
> $$\Psi_F(t) := 4\,\mathfrak{m}_F\,(e^{t/2}+e^{-t/2}-2) - \sum_{n \leq e^t} \frac{c_F(n)}{\sqrt{n}}(t-\log n) + \left[ \log Q + \sum_{j=1}^{r} \lambda_{j} \frac{\Gamma'}{\Gamma}\left(\frac{\lambda_{j}}{2} +\mu_{j}\right) \right] t + \sum_{j=1}^{r} \lambda_{j}^2 \left[\Phi\left(1,2,\frac{\lambda_{j}}{2}+\mu_{j}\right)- e^{-t(\frac{1}{2}+\frac{\mu_{j}}{\lambda_{j}})}\Phi\left(e^{-\frac{t}{\lambda_{j}}},2,\frac{\lambda_{j}}{2}+\mu_{j}\right)\right]$$
> for non-negative $t$ [...] Then $\Psi_F(t)$ is uniquely determined from $F$ and $\Psi_F(t)=-g_F(t)$ holds for all real numbers $t$.

Here `c_F(n)` are the coefficients of `-F'/F(s) = sum c_F(n) n^{-s}` (axiom (S5')), `m_F` the order of the pole at `s=1`, and `Phi` the Lerch transcendent.

What this is and is not, relative to verdict (iii):

- It is a GRH equivalence with eventual one-signedness of a `log(x/n)`-smoothed sum, valid for a class that contains the normalised `L(s + 11/2, Delta)`.
- The smoothed sum runs over **prime powers** with von Mangoldt-type coefficients `c_F(n)`, and the quantity whose sign is tested is the full screw function: the arithmetic sum **plus** an explicit linear term in `t = log x` and a bounded archimedean term. It is not a statement about the sign of a sum of `tau(p)` over primes alone.
- The hypothesis "no real zeros except possibly at `s = 1/2`" must be checked for the `F` in question. **(compiler's remark)** For Delta it follows from positivity of `Delta(iy)` (the completed L-function is a Mellin transform of a positive function, as used in Koyama-Sheth Theorem 3.9); I did not find this specialisation written in Suzuki's paper.
- `t_0` is existential.

---

## 7. Rubinstein-Sarnak and Sarnak's letter

### 7.1 Rubinstein-Sarnak 1994  **[P]** (pages 173-176: abstract, Section 1 with Theorems 1.1 and 1.2 and the GSH statement; pages 193-197: end of Section 4, Section 5 "Generalizations", references. Sections 2-3 and most of Section 4 not read.)

M. Rubinstein, P. Sarnak, "Chebyshev's Bias", Experimental Mathematics 3 (1994), no. 3, 173-197. Read from the copy on M. Rubinstein's web page (page images; the PDF has no usable text layer).

> **Theorem 1.1** (see Section 2.1). *Assume GRH. Then* $E_{q;a_1,\dots,a_r}$ *has a limiting distribution* $\mu_{q;a_1,\dots,a_r}$ *on* $\mathbb{R}^r$, *that is,* $\lim_{X\to\infty}\frac{1}{\log X}\int_2^X f(E_{q;a_1,\dots,a_r}(x))\frac{dx}{x}=\int_{\mathbb{R}^r} f(x)\,d\mu_{q;a_1,\dots,a_r}(x)$ *for all bounded continuous functions* $f$ *on* $\mathbb{R}^r$.

> **Theorem 1.2** (see Section 2.2). *Assume GRH. There are positive constants* $c_1, c_2, c_3$ *and* $c_4$, *depending only on* $q$, *such that* $\mu_{q;a_1,\dots,a_r}(B'_R)\le c_1\exp(-c_2\sqrt R)$, $\mu_{q;a_1,\dots,a_r}(B^{\pm}_R)\ge c_3\exp(-\exp c_4 R)$.

(Transcribed by me from the page images.) Also on p. 174: "Note, however, that assuming only GRH we don't know that $\delta(P_{q;a_1,\dots,a_r})$ exists". Under GRH and GSH: `delta(P_{4;3,1}) = 0.9959...`, `delta(P_{3;2,1}) = 0.9990...`, and for `pi(x)` versus `Li(x)`, `delta(P_1) = 0.00000026...`.

**What it says about tau: nothing.** The paper concerns `pi(x;q,a)` races and `pi(x) - Li(x)`. Section 5 ("Generalizations") discusses prime ideals in ideal classes, Chebotarev-type questions for Galois extensions (with a caution about GSH for Artin L-functions and Armitage's example of a zero at `s = 1/2`), and prime geodesics in homology classes. Tau, modular forms and elliptic curves do not appear in the parts I read, and the reference list contains no modular-form item.

### 7.2 Sarnak's letter to Mazur  **[F]**

P. Sarnak, "Letter to: Barry Mazur on 'Chebyshev's bias' for tau(p)", dated November 2007, with an addendum "Added - January, 2008"; the file is stamped "August 12, 2008". Unpublished letter; read in full (11 pages) from Sarnak's Princeton web page. Equation numbers below are the letter's.

Hypotheses (verbatim):

> (a) $L(s,\pi,\mathrm{sym}^n)$ extends to an entire function in $s$ and satisfies an expected functional equation [...]
> (b) $L(s, \pi, \mathrm{sym}^n)$ satisfies the Riemann Hypothesis so that $\rho_{j,n} = \frac12 + i\gamma_{j,n}$ with $\gamma_{j,n}\in\mathbb{R}$.
> (c) The numbers $\gamma_{n,j} > 0$, $n \ge 1$ and $j \ge 1$ are linearly independent over $\mathbb{Q}$.

(c) is introduced as "the working hypothesis". Under (a) and (b), for smooth `V` on `[0, pi]` orthogonal to `U_0`, with `A_n = <V, U_n>`: `S_V(x) := (log x / sqrt x) sum_{p<=x} V(theta_p)` "has a limiting distribution $\mu_V$, w.r.t. $dx/x$" (18)-(19), with mean (20)
`E(S_V) = sum_{n>=1} A_n ( -2 M_n(1/2) + (-1)^{n+1} )`, where `M_n(1/2)` is the multiplicity of the zero of `L(s, pi, sym^n)` at `1/2`. Under (a), (b), (c): Fourier transform formula (21), and "the measure $\mu_V$ is symmetric about its mean $E(S_V)$, it is smooth and has support all of $(-\infty,\infty)$", with variance (23) `V(S_V) = sum_{n>=1} sum_{gamma_{j,n}>0} |A_n|^2 / (1/4 + gamma_{j,n}^2)`.

Source of the bias (verbatim): "First we examine the second term in (7) which is one source of the 'bias'." That second term is the sum over prime squares `sum_{p^2<=x} log p U_n(2 theta_p)` in the explicit formula for `L'/L(s, pi, sym^n)`; using `<U_m(2 theta), U_0> = (-1)^m` (12) it contributes `(-1)^{n+1} x^{1/2}` (13). The other source is the central zeros `M_n(1/2)` (14). For `n = 1` the prime-square term is governed by `U_1(2 theta) = alpha^2 + beta^2 + 1`, i.e. by the fact that `Sym^2` has no pole at `s = 1`.

For tau, weighted sum (verbatim):

> Applying this to $V(\theta) = U$ [sic], yields a limiting distribution, $\mu_1$ for
> $$S_1(x) := \frac{\log x}{\sqrt{x}}\sum_{p\le x}\lambda_p. \quad (24)$$
> Now $m_1\left(\frac12\right) = 0$ so
> $$E(S_1) = 1. \quad (25)$$
> Thus the weighted sum over $\tau(p)/p^{11/2}$ has a definite bias to being positive.

So for `S_1` the letter gives: existence of the limiting logarithmic distribution and mean 1 under (a)+(b) (only `n = 1`, i.e. `L(s, Delta)`, and the non-vanishing at `s=1` of `Sym^2` enter); symmetry, smoothness, full support and the variance formula under (c) as well. The letter gives no numerical value of the density `delta({x : S_1(x) > 0})` for Delta.

For the sign race `#{p<=x : tau(p)>0} - #{p<=x : tau(p)<0}` (Mazur's Figure 1.6), with `H` the Heaviside function on `[0, pi]`: the letter needs all symmetric powers, extends the working hypothesis to `M_n(1/2) = (1 - epsilon(1/2, pi, sym^n))/2` (33), with root number `+1` for `n = 0,1,2,4,6,7 mod 8` and `-1` for `n = 3,5 mod 8` (34), and finds mean (41)
`E(S) = (2/pi) (1 + 2 int_0^1 x^2/(1+x^4) dx)`, while (42)-(43) give `V(S_N) >> log N`. Verbatim conclusion:

> This indicates that for $S$ in (38) $V(S)=\infty$. That is for your fluctuations in (1.6), while the mean is positive the fluctuations are big enough to dissolve the bias. Put another way $\delta(\{x : S(x) > 0\}) = \frac12$.

Elliptic curves: `S_{1,E}(x) = (log x/sqrt x) sum_{p<=x} a_E(p)/sqrt p` has mean `-2 r(E) + 1` (27) and variance `sum_{gamma_{j,E}>0} 1/(1/4+gamma_{j,E}^2)` (28).

January 2008 addendum (elliptic curves only): for `G(x) = sum_{p<=x} a_E(p)/p`, "$G(x) = r \log\log x + c_4(E) + o(1)$ iff $S(x) = o(\log x)$ (xi)", where `S(x) = sum_{rho != 1/2} x^{i gamma}/(i gamma)`; "It is hard to imagine that the second condition in (xi) is not true". This is the `p^{-1/2}`-extra-weight sum (the Aoki-Koyama weighting) and the condition is of DRH type; it is stated for `E`, not for tau.

Status: a letter with derivations at the level of "it follows as in the analysis in [2]"; predictions rather than theorems with full proofs. The rigorous general framework is Devin's (section 8).

---

## 8. Devin, Chebyshev's bias for analytic L-functions  **[P]** (arXiv v4 TeX: introduction, Sections 1 and 2 in full; the subsection "L-functions of automorphic forms on GL(m)" containing Theorem 3.4; Lemma 5.7 and Corollary 5.8 with its proof. Section 4 (proof of Theorem 2.1), most of Section 5 and the other examples in Section 3 not read.)

L. Devin, "Chebyshev's bias for analytic L-functions", arXiv:1706.06394 (v1 20 Jun 2017, v4 26 Apr 2018). Math. Proc. Cambridge Philos. Soc. 169 (2020), no. 1, 103-140 (published online 22 March 2019). Theorem numbers below are those of the arXiv v4 PDF; the published version was not compared.

Setting: "analytic L-functions" (Definition 1.1): Euler product of degree `d` with `|alpha_j(p)| = 1` for `p` not dividing the conductor (Ramanujan-Petersson built in), functional equation, and a second-moment hypothesis. Example 1(3): "Modular $L$-functions of degree $2$ are analytic $L$-functions. It is a consequence of results of Deligne and Serre". So the normalised `L(s + 11/2, Delta)` is covered unconditionally by the framework. The race studied is the **unweighted** one in the normalised coefficients, `sum_{p<=x} lambda_f(p)`, i.e. `sum_{p<=x} tau(p) p^{-11/2}` for Delta; not the `p^{-6}` sum.

> **Theorem 2.1.** [...] Define
> $$E_{\mathcal{S}}(x)=\frac{\log x}{x^{\beta_{\mathcal{S},0}}}\left(\sum_{p\leq x}\sum_{f\in\mathcal{S}}a_{f}\lambda_{f}(p) + \sum_{f\in\mathcal{S}}a_{f}m(L(f,\cdot),1) \Li(x) \right).$$
> The function $E_{\mathcal{S}}(x)$ admits a limiting logarithmic distribution $\mu_{\mathcal{S}}$. There exists a positive constant $C$ (depending on $\mathcal{S}$) such that one has $\mu_{\mathcal{S}}(\mathbf{R}-[-A,A]) \ll \exp(-C\sqrt{A}).$ Moreover let $X_{\mathcal{S}}$ be a random variable of law $\mu_{\mathcal{S}}$, then the expected value of $X_{\mathcal{S}}$ is
> $$\mathbb{E}(X_{\mathcal{S}}) = m_{\mathcal{S}} := \sum_{f\in\mathcal{S}}a_{f}\left(m(L(f^{(2)},\cdot),1)\delta_{\beta_{\mathcal{S},0}=1/2} -\beta_{\mathcal{S},0}^{-1}m(L(f,\cdot),\beta_{\mathcal{S},0})\right),$$
> and its variance is
> $$\Var(X_{\mathcal{S}}) = 2\sum_{\gamma \in \mathcal{Z}_{\mathcal{S}}^{*}}\frac{\lvert M(\gamma)\rvert^{2}}{\beta_{\mathcal{S},0}^{2}+\gamma^{2}}$$

This is unconditional (no GRH, no LI); `beta_{S,0}` is the supremum of real parts of zeros.

> **Theorem 2.5.** Suppose that the set $\mathcal{S}$ satisfies the conditions of Theorem~\ref{Th_DistLim}. (1) Suppose the Riemann Hypothesis is satisfied for every $L(f,s)$, $f\in\mathcal{S}$ (\textit{i.e.} $\beta_{\mathcal{S},0}=\frac{1}{2}$). Suppose also that for every $f\in \mathcal{S}$, one has $\re(a_{f})\geq 0$, and that there exists $f\in\mathcal{S}$ such that $\re(a_{f})> 0$. Then there exists a constant $c$ depending on $\mathcal{S}$ such that
> $$\mu_{\mathcal{S}}(\mathbf{R}-[-A,A]) \gg \exp(-\exp(cA)).$$
> In particular $0< \underline{\delta}(\mathcal{S}) \leq \overline{\delta}(\mathcal{S}) <1$. (2) Suppose $\beta_{\mathcal{S},0}>\frac{1}{2}$, and the sum $\sum_{\gamma \in \mathcal{Z}_{\mathcal{S}}}\frac{1}{\lvert \beta_{\mathcal{S},0} + i\gamma\rvert}$ converges. Then $\mu_{\mathcal{S}}$ has compact support.

where (Definition 1.3) `overline-delta(S)` and `underline-delta(S)` are the upper and lower logarithmic densities of `{x : S(x) >= 0}`. The remark after the theorem says the hypothesis on `a_f` "is satisfied in the case the set $\mathcal{S}$ is a singleton".

> **Corollary 5.8.** In the case $m_{\mathcal{S}}<0$ one has: $\overline{\delta}(\mathcal{S}) \leq \frac{2}{m_{\mathcal{S}}^{2}}\sum_{\gamma \in \mathcal{Z}_{\mathcal{S}}^{*}}\frac{\lvert M(\gamma)\rvert^{2}}{(\beta_{\mathcal{S},0}^{2}+\gamma^{2})}.$ In the case $m_{\mathcal{S}}>0$ one has: $\underline{\delta}(\mathcal{S}) \geq 1- \frac{2}{m_{\mathcal{S}}^{2}}\sum_{\gamma \in \mathcal{Z}_{\mathcal{S}}^{*}}\frac{\lvert M(\gamma)\rvert^{2}}{(\beta_{\mathcal{S},0}^{2}+\gamma^{2})}.$

(Proved from Theorem 2.1 and Chebyshev's inequality, "inspired from \cite[Lem. 2.7]{FioEC}"; no hypothesis on linear independence.)

> **Theorem 3.4.** Let $L(\pi,s)$ be a real $L$-function associated to an irreducible unitary cuspidal automorphic representation of $GL(m)$ with $m\geq 2$. Suppose the Ramanujan--Petersson conjecture holds for $L(\pi,s)$. Then [...] the function $E_{\pi}(x) = \frac{\log x}{x^{\beta_{\pi,0}}}\sum_{p\leq x} \lambda_{\pi}(p)$ has a limiting logarithmic distribution $\mu_{\pi}$. Moreover under GRH for $L(\pi,s)$, the mean of $\mu_{\pi}$ is $m_{\pi} = \pm 1 -2m(L(\pi,\cdot),1/2)\neq 0$.

followed by "Morally, under GRH, since the mean value is non zero, we expect that the prime number race associated to such an $L$-function has always a bias." and "In the case $m=2$, the Ramanujan--Petersson conjecture has been proved by works of Deligne and Deligne--Serre".

Existence of the density (Theorem 2.2) needs a partial linear-independence hypothesis ("there exists $N\geq 1$ and $\lambda_{1},\ldots,\lambda_{N} \in \mathcal{Z}_{\mathcal{S}}^{*}$ such that $\langle\lambda_{1},\ldots,\lambda_{N}\rangle_{\mathbf{Q}} \cap \langle \mathcal{Z}_{\mathcal{S}}^{*} - \lbrace \lambda_{1},\ldots,\lambda_{N}  \rbrace \rangle_{\mathbf{Q}} = \lbrace 0 \rbrace$"); Corollary 2.4 then gives `(delta(S) - 1/2) m_S >= 0` under that plus either a parity-type independence condition or the concentration condition `(2/m_S^2) sum |M(gamma)|^2/(beta^2+gamma^2) < 1/2`.

**Answers to the questions asked.**

- Existence of the limiting logarithmic distribution for `sum_{p<=x} lambda_f(p)` for modular forms: yes, and unconditionally (Theorem 2.1 with Example 1(3)), not merely under GRH.
- Positive lower densities under GRH without LI: yes, Theorem 2.5(1) gives `0 < underline-delta <= overline-delta < 1` for a singleton `S`. So under GRH alone both `{S(x) >= 0}` and its complement have positive lower logarithmic density, for the unweighted normalised sum.
- A quantitative lower bound under GRH without LI: Corollary 5.8 (Chebyshev's inequality), in terms of `sum_gamma |M(gamma)|^2/(1/4+gamma^2)`.
- **Tau specifically: my recollection of a worked tau example was wrong.** In the arXiv v4 text the only occurrence of Ramanujan's tau is in the description of Sarnak's letter: "In \cite{SarnakLetter}, Sarnak presents a method to deal with this question in the case $\mathcal{S}$ is a singleton. Precisely he considers the cases of the $L$-function associated to Ramanujan's $\tau$ function and of $L$-functions associated with an elliptic curve over $\mathbf{Q}$." (plus the title of the letter in the bibliography). The worked examples in Section 3 are: `pi(x)` vs `Li`, races mod `q`, automorphic `GL(m)` (Theorem 3.4), elliptic curves (Proposition 3.5), sums of two squares / Gaussian primes (Theorems 3.6, 3.7), correlations of two elliptic curves (Theorem 3.8), and `J_0(q)` (Theorem 3.11). No mean, variance or density is computed for Delta. **(compiler's remark)** Specialising Theorem 3.4 / Theorem 2.1 to Delta under GRH gives mean `+1` (second-moment L-function has a simple zero at `s=1`, `L(1/2) != 0`), which agrees with Sarnak's (25); the variance sum over the zeros of `L(s, Delta)` is not evaluated in any source I read.

Cited by Devin and not read by me **[N]**: Akbary-Ng-Shahabi, "Limiting distributions of the classical error terms of prime number theory" (cited for `psi(f,x)/sqrt x` under GRH without Ramanujan-Petersson, their Cor. 1.5); Kaczorowski-Ramare, Acta Arith. 106 (2003) 277-297 (Devin: they "prove the existence of the logarithmic density of a comparable set in the general setting of the Selberg class, under the Riemann Hypothesis only"); Fiorilli's elliptic-curve paper [FioEC].

---

## 9. Searches (item 9)

### 9.1 arXiv API queries (export.arxiv.org/api/query), run 2026-11-24

Results are complete lists unless a count is given. The API's phrase matching is on title/abstract/metadata, so papers whose abstracts avoid the search words are missed.

| Query (abbreviated) | Hits | Relevant hits |
|---|---|---|
| `"Chebyshev" AND "bias" AND ("Ramanujan" OR "tau function" OR "modular forms" OR "cusp form" OR "Hecke eigenvalues")` | 4 | 2509.04187, 2405.01512, 2203.12791 (4th: Moree 2001, unrelated) |
| `("prime number race(s)" OR "prime races") AND ("Fourier coefficients" OR "cusp form" OR "modular form" OR "automorphic")` | 0 | - |
| `"Deep Riemann Hypothesis"` | 10 | 2611.20592 (Aoki), 2607.28931, 2607.23150, 2307.05958 (Okumura, Fermat curves), 2206.05445 (Kaneko-Koyama, function fields), 2206.02612 (Kaneko-Koyama-Kurokawa), 2203.12791, 2203.12266, 1902.04203 (Kaneko), 1210.1216 (Kimura-Koyama-Kurokawa) |
| `("Chebyshev's bias" OR "Chebyshev bias" OR "prime number race(s)")`, 2023-2026 | 16 | as above plus 2608.06435, 2603.21005, 2603.20093, 2512.23302, 2505.02968, 2412.19728, 2412.14792, 2404.06804, 2401.04000, 2308.13959, 2302.13665; none on tau or modular forms beyond those already listed |
| `Ramanujan AND tau AND primes AND Riemann` (abstract), math.NT, 2022-2026 | 0 | - |
| `"Hecke eigenvalues" AND primes AND "Riemann Hypothesis"` (abstract), 2022-2026 | 0 | - |
| `("Hecke eigenvalues" OR "Fourier coefficients") AND explicit AND primes AND (GRH)`, 2022-2026 | 0 | - |
| `"Fourier coefficients" AND cusp AND primes AND (oscillation OR "sign changes" OR Landau OR omega)`, 2022-2026 | 2 | neither about sums over primes of `tau(p)` |
| `("equivalent to the (generalized) Riemann hypothesis" OR "equivalent to GRH") AND (modular OR automorphic OR Ramanujan OR "cusp form" OR "Selberg class")` | 1 | unrelated |
| `"logarithmic measure" AND ("Euler product" OR bias)` | 3 | 2512.23302, 2405.01512, 2312.05236 |
| `(Cesaro OR "Riesz mean(s)" OR smoothed) AND ("Hecke eigenvalues" OR "Ramanujan tau" OR "Fourier coefficients of cusp forms") AND primes` | 4 | none relevant (1305.5283 Rouse-Thorner uses a smooth cutoff for Sato-Tate counting) |
| `("Ramanujan" AND "tau") AND ("logarithmic density" OR "limiting distribution" OR "sign changes" OR oscillat*) AND prime` | 0 | - |
| `"symmetric square" AND zeros AND (rigorous OR verification OR "interval arithmetic")` | 0 | - |
| `"symmetric power" AND "L-functions" AND (compute OR computation OR numerical)` | 14 | 0604095 (Martin-Watkins, elliptic curves; not read) |
| title searches 2023-2026: `ti:bias`, `ti:race(s)`, `ti:Chebyshev`, `ti:"tau function"`/`"Ramanujan tau"` (math.NT) | 10 / 6 / 24 / 6 | nothing new on tau prime sums |
| author searches: Sheth, Hayani, Suzuki | - | found 2209.12832, 2611.12152 |

### 9.2 Web searches

Four web searches (terms: Chebyshev's bias + Ramanujan tau + GRH + sum over primes; equivalence of RH with positivity / constant sign of smoothed prime sums for tau, Hecke eigenvalues or automorphic L-functions; zeros and rigorous verification for the symmetric square of Delta; rigorous verification for modular-form L-functions). They returned the papers already listed, plus Yoshida's paper (9.3) and Han's (5.3). No hit for any explicit, all-`x >= x_0` result on sums over primes of `tau(p)` or Hecke eigenvalues under GRH.

### 9.3 Explicit results at primes for tau / Hecke eigenvalues found  **[A]**

Nearest items with explicit constants (abstracts only):

- J. Rouse, J. Thorner, "The Explicit Sato-Tate Conjecture and Densities Pertaining to Lehmer-Type Questions", arXiv:1305.5283, Trans. Amer. Math. Soc. 369 (2017) 3575-3604. Abstract: assuming analytic properties of the symmetric power L-functions "(all of which follow from Langlands functoriality and the Generalized Riemann Hypothesis), we prove that if $x$ is sufficiently large, then $\left|\#\{p\leq x:\theta_p\in I\} -\mu_{ST}(I)\int_2^x\frac{dt}{\log t}\right|\ll\frac{x^{3/4}\log(N k x)}{\log x}$ with an implied constant of $3.34$", and for tau "$\lim_{x\to\infty}\frac{\#\{n\leq x:\tau(n)\neq0\}}{x}>1-1.54\times10^{-13}$". Explicit constant, but GRH for all symmetric powers, "x sufficiently large", and a counting statement, not a sign statement for `sum tau(p) p^{-s}`.
- A. Hoey, J. Iskander, S. Jin, F. Trejos Suarez, arXiv:2108.03520: unconditional explicit Sato-Tate error `58.1 log((k-1)N log x)/sqrt(log x)` for `x >= 3`. Not a prime-sum sign statement.
- J. Thorner, arXiv:1407.2656 (Arch. Math. 103 (2014)): effective Sato-Tate error term assuming automorphy of symmetric powers.

I found no paper giving numerical constants for `psi(x, Delta)`-type sums or for `sum_{p<=x} tau(p) p^{-s}` under GRH.

### 9.4 Zeros of `L(s, Sym^2 Delta)` and rigorous verification

- H. Yoshida, "On calculations of zeros of various L-functions", arXiv:math/9411213, J. Math. Kyoto Univ. 35 (1995), no. 4 **[P]** (introduction and Section 6 read). First sentence: "As we have shown several years ago [Y2], zeros of $L(s,\Delta)$ and $L^{(2)}(s,\Delta)$ can be calculated quite efficiently by a certain experimental method." The symmetric-square zeros themselves are in [Y2] = H. Yoshida, "On calculations of zeros of L-functions related with Ramanujan's discriminant function on the critical line", J. Ramanujan Math. Soc. 3 (1988), 87-95 **[N]** (not obtained; the 1995 paper treats only the 3rd and 4th symmetric powers of Delta, "Since the cases $r = 1, 2$ are discussed in [Y2]"). The method is explicitly non-rigorous: "The most serious defect of our method of calculation is that we do not have rigorous controle of error estimates." So these are numerical zeros on the critical line, not a verification of RH up to a height. I did not obtain the numerical values.
- A. Best, J. Bober, A. Booker, E. Costa, J. Cremona, M. Derickx, M. Lee, D. Lowry-Duda, D. Roe, A. Sutherland, J. Voight, "Computing classical modular forms", arXiv:2002.04717 **[P]** (only the L-function passages were read, located by text search). For the **degree-2** L-functions of the newforms in the LMFDB: computations use "rigorous error bounds and interval arithmetic", and "As a consistency check for our $L$-function computations, after computing a provisional list of all non-trivial zeros on the critical line up to a chosen height bound $b$ we confirmed that no zeros are missing, in other words, that the Riemann Hypothesis holds for each $L$-function up to height $b$. We use the method described in \cite{buthe2015} based on the Weil--Barner explicit formula." The value of `b` is not in the passages I read. The underlying general method is deferred to "Jonathan W. Bober, Andrew R. Booker, Edgar Costa, Min Lee, David J. Platt, and Andrew Sutherland, *Computing motivic L-functions*, in preparation"; I did not find that paper on arXiv. This covers `L(s, Delta)` (newform 1.12.a.a), not `Sym^2 Delta`.
- A. Booker, "Artin's conjecture, Turing's method, and the Riemann hypothesis", Experiment. Math. 15 (2006) **[A]** (search snippet only): a rigorous algorithm for general L-functions on the critical line and a generalisation of Turing's method, applied to Dedekind zeta and Artin L-functions. I did not check whether any symmetric square is treated there.
- LMFDB, what was actually obtained:
  - `https://www.lmfdb.org/L/3/1/1.1/r0e22/` and `https://www.lmfdb.org/L/SymmetricPower/2/ModularForm/GL2/Q/holomorphic/1/12/a/a/`: reCAPTCHA page via curl; **HTTP 404** via the web-fetch tool.
  - Search `https://www.lmfdb.org/L/?degree=3&conductor=1&motivic_weight=22` and `...&algebraic=rational`, via the web-fetch tool: reported "No matches" (0 results) for both. Caveat: that tool returns a model-written summary of the page, and shortly afterwards the site began returning a browser-check page, so I could not re-verify. The one successful JSON API call (first 100 of the degree-3, conductor-1 records) contained only `ModularForm/GL3/Q/Maass` origins.
  - The newform page `1.12.a.a` lists as related objects only the newspaces and "L-function 1.12.a.a"; no symmetric-square link (same caveat).
  - Net: **I could not find a `Sym^2 Delta` L-function page in the LMFDB and obtained no zeros of `L(s, Sym^2 Delta)` from it.** I cannot rule out that it exists under a label I did not try.
  - For `L(s, Delta)` itself (label `2-1-1.1-c11-0-0`), via the web-fetch tool: `Lambda(s) = Gamma_C(s + 11/2) L(s) = Lambda(1 - s)`; `L(1/2) ~ 0.7921228386`; first ten zeros (imaginary parts) 9.222379399921102522243767192743, 13.90754986139213440644668132877, 17.44277697823447331355152513713, 19.65651314195496100012728175632, 22.33610363720986727568267445924, 25.27463654811236535674532419313, 26.80439115835040303257574923358, 28.83168262418687544502196191298, 31.17820949836025906449218889077, 32.77487538223120744183045567331. (Relayed by the summarising tool; not checked against raw page data.) The LMFDB knowl "Reliability of L-function data" (`rcs.rigor.lfunction`) was reported as containing only "TBA".

---

## 10. Do results of the exact types (a), (b), (c) exist?

**(a) "Under GRH for L(s, Delta), `sum_{p<=x} tau(p) p^{-11/2} (1 - p/x) > 0` for every `x > x_0` (explicit)".** Not found, neither with an explicit `x_0` nor in qualitative form, for tau or for any modular form. Nearest:

- Hayani Theorem 1.4 (section 4): the same Cesaro smoothing, GRH alone, every `x`, but for Dirichlet characters and the weight `p^{-1/2}`, with a non-explicit `O`-term.
- Suzuki Theorems 3, 4, 8, 9 (section 6.1): `log(x/n)` smoothing, Dirichlet characters and zeta, existential `x_0`.
- Johnston (zeta): the only item in this file with a fully explicit all-`x` range (`x > 2`).
- Note the normalisation in (a): weight `p^{-11/2}` with a `(1 - p/x)` kernel is a Cesaro mean of Devin's/Sarnak's unweighted `S(x)`, not of the Aoki-Koyama `p^{-6}` sum. I found no source discussing a Cesaro mean of `sum_{p<=x} lambda_f(p)` for a cusp form.

**(b) "GRH for L(s, Delta) is equivalent to the eventual one-signedness of a smoothed sum of tau(p) over primes".** Not found as stated. Nearest:

- Suzuki, IJNT 21 (2025), Theorems 1.1 and 4.1 (section 6.2): an equivalence for the extended Selberg class, hence applicable in principle to `L(s + 11/2, Delta)`, but for the full screw function (prime-power sum with weight `(t - log n)/sqrt n`, plus a linear term and a bounded archimedean term), and without Delta or tau being mentioned.
- Suzuki, Ramanujan J. (2025), Theorems 4 and 9: equivalence with a sum over primes only, for real Dirichlet characters.
- The commented-out converse in Sheth's TeX source (section 3): unpublished.
- Aoki-Koyama / Koyama-Kurokawa / Aoki 2611.20592: equivalences with DRH, not GRH (as described in Hayani's introduction and Aoki's abstract; I did not read those equivalence proofs).

**(c) "Lower bounds for the logarithmic density of `{x : sum_{p<=x} tau(p) p^{-11/2} > 0}` under GRH without LI".** No bound specific to tau found. General results that apply:

- Devin Theorem 2.5(1): under RH for `L(s,f)`, `0 < underline-delta <= overline-delta < 1` (qualitative).
- Devin Corollary 5.8: `underline-delta >= 1 - (2/m^2) sum_gamma |M(gamma)|^2/(1/4 + gamma^2)` when `m > 0` (no LI needed; the formula is stated without GRH, with `beta_0` in place of `1/2`).
- No source evaluates the right-hand side for Delta. Sarnak's letter gives mean 1 and the variance formula under GRH (+LI for the Fourier transform) but no number.

---

## Verdicts

**(i) Is "`sum_{p<=x} tau(p)/p^6 > 0` on a set of natural (or logarithmic) density 1 under GRH alone" already in the literature?**

In substance, yes. Sheth, Math. Proc. Camb. Phil. Soc. 179 (2025) 331-349, Corollary B = Corollary 4.3 (arXiv:2405.01512v1: Corollary 4.4): assuming only RH for `L(s, Delta)`, `sum_{p<=x} tau(p)/p^6 = (1/2) log log x + c + o(1)` for all `x` outside a set of finite logarithmic measure. Positivity off a set of finite logarithmic measure is immediate, and natural density 1 follows from the elementary remark written out in Hayani, arXiv:2512.23302, display (1.5). Koyama-Sheth, arXiv:2509.04187, repeat the GRH-only remark for all newforms after their Theorem 1.5.

Caveats, all checked against the texts: (1) the literal sentence "natural density 1 under GRH" for tau was not found; Koyama-Kurokawa's Corollary 1 states natural density 1 under **DRH**, Hayani's Corollary 1.2 states it under GRH for the mod 4 race only. (2) Sheth's `o(1)` has no rate and the exceptional set is not described; the rate `(log log x)^{3+eps}/log x` exists only in Hayani's Dirichlet setting. (3) Sheth's key steps are proved by citation (Gallagher/Qu for the exceptional set, Sheth 2023 for the Euler-product formula), which I did not read.

**(ii) Are any all-`x` / explicit-threshold statements for tau known?**

None found. For tau the only all-`x` asymptotic is Koyama-Kurokawa's Theorem 2, under DRH and with no effectivity. Under GRH alone, all-`x` statements exist only outside the modular-form setting: Hayani Theorem 1.4 (Cesaro mean, Dirichlet characters, non-explicit constant), Suzuki's theorems (Dirichlet characters and zeta, existential `x_0`), Johnston (zeta, `x > 2`), and the special Artin-type extreme-bias races in Hayani 2607.23150. No source gives a numerical `x_0` for any prime sum of `tau(p)`.

**(iii) Is the equivalence "GRH(Delta) iff eventual sign of a smoothed prime sum of tau" stated anywhere?**

Not for Delta or tau, and not for a sum over primes of a cusp form's coefficients. The closest published statement is Suzuki's Theorem 1.1 (with the zero-free formula of Theorem 4.1) in IJNT 21 (2025) 1815-1823 = arXiv:2209.12832: for any `F` in the extended Selberg class with no real zeros other than possibly `s = 1/2`, GRH for `F` iff the real part of the negated screw function is nonnegative for all `t >= t_0`. That function is a `log(x/n)`-smoothed prime-power sum plus explicit archimedean terms, and the paper names "L-functions of holomorphic cusp forms" as members of the class but does not specialise. Pure prime-sum equivalences (Suzuki, Ramanujan J. 2025, Theorems 4 and 9) are stated for real Dirichlet characters only. A converse "bounded central Euler product implies RH" appears only in commented-out draft text in Sheth's arXiv source.

**(iv) Is a rigorous RH verification for `Sym^2 Delta` known?**

Not found. Zeros of `L(s, Sym^2 Delta)` were computed by Yoshida (J. Ramanujan Math. Soc. 3 (1988), not obtained; described in his 1995 paper) by a method he himself calls non-rigorous. Rigorous interval-arithmetic computations with RH checked up to a height exist for the degree-2 modular-form L-functions in the LMFDB, including `L(s, Delta)` (Best et al., arXiv:2002.04717; height bound not stated in the passages read). I could not locate a `Sym^2 Delta` page in the LMFDB (old URLs 404, a degree-3/conductor-1/motivic-weight-22 search reported no matches, then the site blocked further access), and found no paper reporting a rigorous verification for this degree-3 L-function. This is an absence of evidence from a limited search, not a proof that no such computation exists; in particular I did not read Booker 2006 and could not find the announced "Computing motivic L-functions".

---

## Appendix: things that were not read and could change a verdict

- Published version of Koyama-Kurokawa (Proc. Japan Acad.), in case it differs from arXiv v2.
- Gallagher 1980, Qu 2007, Sheth arXiv:2312.05236, Conrad 2005, Goldfeld 1982 (inputs to Sheth's theorems).
- Kaneko-Koyama-Kurokawa, arXiv:2206.02612, Section "Applications to the Chebyshev Bias" (only the abstract and a keyword search were done); their announced "Chebyshev's bias for GL_n".
- Kaczorowski-Ramare 2003 and Akbary-Ng-Shahabi 2014 (GRH-only limiting distributions / densities in the Selberg class).
- Rubinstein-Sarnak Sections 2-4; Devin Sections 3 (other examples), 4, 5.
- Proofs in Suzuki's two papers; Section 5 of arXiv:2209.12832 ("generalizations of other results").
- Yoshida 1988; Booker 2006; any LMFDB data for symmetric squares.
- The MathSciNet / zbMATH citation lists of Koyama-Kurokawa and Sheth were not consulted (no access attempted), so papers citing them that are not on arXiv would have been missed.
