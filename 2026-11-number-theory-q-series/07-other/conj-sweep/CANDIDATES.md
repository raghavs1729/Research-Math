# conj-sweep: candidates

Harvest (26 Nov 2026): arXiv API (`export.arxiv.org/api/query`), 15 queries over math.NT / math.CO /
math.CA with submission dates 2024-01-01 … 2026-11-26 (phrases such as "we conjecture", "numerical
evidence", "several conjectures", topic words: supercongruence, continued fraction, log-concave,
partition congruences, Hankel, OEIS, positivity, q-series), ≈ 2900 abstracts; 358 abstracts announce
conjectures after filtering out graph theory / algebraic geometry / analytic number theory of
L-functions; for 65 of them the TeX source was read and the `conjecture` environments extracted
(verbatim below, truncated at 520 characters; sources kept outside the repository).
Excluded on instruction: Nahm sums / Kanade–Russell identities, Ramanujan Machine constants, Zudilin
L-values, Ramanujan tau, Sun's series, Erdős problems.

## Status of each candidate

| arXiv id | status | where |
|---|---|---|
| 2603.08510 | all four conjectural congruences PROVED | R02 |
| 2611.05623 | Conj. 3, 4 PROVED; Conj. 1, 2 survived | R03 |
| 2503.06228 | Conj. 1 PROVED; Conj. 2 (infinite family) not attempted | R03 |
| 2606.16828 | main conjecture, meld and stack conjectures PROVED; Grundy-value conjectures and `conj:hookswap` not settled (no finite Nerode index found) | R07 |
| 2609.16843 | Conj. `conj:except` DISPROVED as printed, corrected and PROVED; Conj. 1, 3 (general periodicity) not attempted | R06 |
| 2605.12822 | Conj. 4: necessity DISPROVED, sufficiency survived; Conj. 1 (q-Fibonomial unimodality) not attempted | R01 |
| 2502.21050 | Conj. 1, 2 DISPROVED as printed (signs); Conj. 3 not tested | R05 |
| 2601.16626 | "if" PROVED, "only if" open | R04 |
| 2506.03620 | both survived | RESULTS T1 |
| 2609.05358 | `conj:period` survived; Conj. 1, 3 (admissible angles) not tested | RESULTS T2 |
| 2404.18517 | Conj. 1–3 survived to n = 700 | RESULTS T3 |
| 2607.16608 | all closing conjectures PROVED | R08 |
| 2611.05196 | survived (p, α) up to (3,4), (5,3), (7,2) | RESULTS T4 |
| 2501.14640 | two KING conjectures survived | RESULTS T5 |
| 2610.06925 | `g ≤ 2h_min` survived; period conjecture not tested | RESULTS T6 |
| 2506.04991 | `conj:heavychoppedrect` DISPROVED; `conj:shallowstair` survived | R09 |
| 2601.18715 | duality conjecture survived | RESULTS T7 |
| 2401.01385 | integral identity confirmed numerically to 68 digits; not proved | RESULTS T8 |
| 2609.01562 | `conj:496` PROVED | R10 |
| 2503.17571 | Conj. 4.1 / 2.6 survived L ≤ 160 | RESULTS T9 |
| 2406.10987 | survived k ≤ 400, n ≤ 4000 | RESULTS T10 |
| 2502.06032 | dropped: authors (Krattenthaler et al.) verified `n ≤ 400` with days of computing; no cheap extension | — |
| 2607.08366 | dropped: certified by the author for `n ≤ 7`, solver search to `n ≤ 13` | — |
| 2609.11813 / 2609.05302 | dropped: the `cφ_18` congruences are the subject of the first paper (already treated there) | — |
| 2409.02789, 2402.15886, 2406.13790, 2609.14295, 2509.10023, 2412.10214, 2602.24239, 2602.22129 | dropped: positivity / asymptotic / algebraic statements without a cheap decisive test | — |
| 2408.13346, 2604.17424, 2406.10987, 2503.17571, 2401.09497, 2508.00946, 2511.19744, 2604.17847, 2609.01562, 2404.18237, 2402.02199, 2512.24221, 2401.01385, 2510.22385, 2610.36132, 2608.08703, 2604.14639, 2609.11072, 2611.28840, 2509.08801, 2508.00074, 2404.15716, 2512.02324, 2512.04775, 2406.02971 | not attempted (budget) | — |

How openness was checked: for the settled items, see `NOVELTY-1.md` (arXiv version history, Semantic
Scholar citations, web search). For items not attempted no search was made.

## Verbatim statements

### arXiv:2605.12822 — Unimodality of $q$-Fibonomial coefficients for small cases
Brendan B. Connelly, Ezekiel Ito, Thomas C. Martinez, Olha Shevchenko, Kacey Yang (2026-05-12). 2 conjecture environment(s) in the TeX source.

* **C1.** `[{\cite[Conjecture 2.5]{kustner}}]\label{conj:mbynunimodal} The polynomials $\qfibonom{m+n}{n}$ are unimodal.`
* **C2.** `\label{conj:unimodalitylemma} Let $\pow \geq2, k \geq 1$ and let $a_1,\ldots, a_k,b$ be positive integers. Then if $\pow \mid a_i$ for any $1 \leq i \leq k$ or \[ b \leq 1+\sum_{i=1}^{k}\left\lfloor\frac{a_i}{\pow}\right\rfloor, \] the polynomial $[a_1]_q\dots[a_k]_q[b]_{q^{\pow}}$ is unimodal. Moreover, if $k \leq 3$ or $\pow \leq 3$, this condition is also necessary.`

### arXiv:2603.08510 — New Ramanujan-type congruences for overpartitions modulo $11$ and $13$
XuanLing Wei (2026-03-09). 2 conjecture environment(s) in the TeX source.

* **C1.** `\begin{equation*} \overline{p}(2^4 n) \equiv 0 \pmod{7}, \quad \text{where } n \equiv 3 \pmod{8} \text{ and } \left(\frac{n}{7}\right)=1. \end{equation*} Equivalently, $\overline{p}(2^4(56n+k)) \equiv 0 \pmod{7}$ for all $n\geq 0$ and $k\in\{11,43,51\}$.`
* **C2.** `\begin{equation*} \overline{p}(17^2\cdot 4n) \equiv 0 \pmod{19}, \quad \text{where } n \equiv 3 \pmod{8},\ \left(\frac{n}{17}\right)=-1,\ \text{and }\left(\frac{n}{19}\right)=1. \end{equation*} It involves $72$ residue classes modulo $2584$; we do not list all of them here.`

### arXiv:2611.05623 — Partitions with designated summands: $k$-tuples and their congruences
Daniel Yaqubi, Madjid Mirzavziri (2026-11-03). 4 conjecture environment(s) in the TeX source.

* **C1.** `\label{conj:A2adic} For every $k\ge1$, the largest power of $2$ dividing all $A_k(4n+3)$ is $2^{v_2(A_k(3))}$ if $k\not\equiv2\pmod{16}$, and exactly $2^3$ if $k\equiv2\pmod{16}$.`
* **C2.** `\label{conj:B2adic} For every $k\ge1$, the largest power of $2$ dividing all $B_k(4n+3)$ is $2^{v_2(B_k(3))}$, except when $k\equiv1\pmod 8$ (where it is exactly $2^2$) or $k\equiv6\pmod 8$ (where it is exactly $2^3$, as established in Theorem~\ref{thm:pdo4exact}).`
* **C3.** `\label{conj:A6} For all $n\ge0$, $A_6(12n+10)\equiv0\pmod{3^5}$, and $3^5$ is the exact power of $3$ dividing all of them.`
* **C4.** `\label{conj:deep} \begin{enumerate} \item[\rm (a)] $A_3(12n+7) \equiv A_3(12n+11) \equiv 0 \pmod{3^4}$. \item[\rm (b)] $A_4(6n+5) \equiv 0 \pmod{3^3}$ and $A_7(12n+11) \equiv 0 \pmod{3^3}$. \item[\rm (c)] $A_8(18n+15) \equiv 0 \pmod{3^5}$. \end{enumerate}`

### arXiv:2503.06228 — Congruences and density results for partitions into distinct even parts
Hemjyoti Nath, Abhishek Sarma (2025-03-08). 2 conjecture environment(s) in the TeX source.

* **C1.** `\label{c1.0} For every $n\geq0$, we have \begin{align*} ped(225n+43) & \equiv 0 \pmod{192},\\ ped(225n+88) & \equiv 0 \pmod{192},\\ ped(225n+133) & \equiv 0 \pmod{192},\\ ped(225n+223) & \equiv 0 \pmod{192}. \end{align*}`
* **C2.** `Let $p\geq 5$ be a prime with $\left( \frac{-2}{p} \right)_{L}=-1$. Also, let $t$ be a positive integer with $(t,6) = 1$ and $p\mid t$. Then for all $n\geq 0$ and $1\leq j \leq p-1$, we have \begin{equation*} ped\left( 9\cdot t^2n + \frac{9\cdot t^2j}{p} + \frac{57\cdot t^2-1}{8} \right) \equiv 0 \pmod{192}. \end{equation*}`

### arXiv:2601.16626 — On generalized eigenvalues of MAX matrices to MIN matrices and of LCM matrices to GCD matrices
Jorma K. Merikoski, Pentti Haukkanen, Antonio Sasaki, Timo Tossavainen (2026-01-23). 1 conjecture environment(s) in the TeX source.

* **C1.** `\label{conj} Let $T=\{1,\dots,n\}$, $n>3$. Then $-1$ is a g-eigenvalue of ${\bf L}_T$ to~${\bf G}_T$ if and only if the binary representation of $n$ begins with $10$.`

### arXiv:2506.03620 — Primes of the form $ax+by$
Yong-Gao Chen, Hui Zhu (2025-06-04). 2 conjecture environment(s) in the TeX source.

* **C1.** `For $3\le a<b$ with $\gcd (a,b)=1$, we have $$\pi (a, b) \ge \frac{13}{66} \pi (s(a,b))$$ and the equality holds if and only if $a=3$ and $b=166$.`
* **C2.** `For $3\le a<b$ with $\gcd (a,b)=1$, we have $$\pi (a, b) \le \frac{1}{2} \pi (s(a,b))$$ and the equality holds if and only if $a=3$ and $b=5$.`

### arXiv:2606.16828 — CRIM: A Natural Game on Integer Partitions
Ina Bašić, Eric Gottlieb, Matjaž Krnc (2026-06-15). 9 conjecture environment(s) in the TeX source.

* **C1.** `\label{conj:main} Every $\p$-position is of even rank.`
* **C2.** `Let $r$ and $k$ be integers with $0 \leq k < r$. Then \[ \sg(R_{r, r}^k)= \begin{cases} 0& r \text{ is even or }k<r-1,\\ 1& r\in \{3,5\}, \, r \text{ is odd, and } k = r-1,\\ 2& \text{otherwise.} \end{cases} \]`
* **C3.** `For $r \geq 7$ we have \[\sg(R_{r,r-1}^k)= \begin{cases} 3& k=r-2 \text{ and } r \text{ odd},\\ 1& \text{otherwise.} \end{cases}\] Also, $\sg(R^1_{3, 2}) = \sg(R^3_{5, 4}) = 3$ and $\sg(R^2_{4, 3}) = 2$ and $\sg(R^4_{6, 5} )= 5$. All other rectairs have Sprague-Grundy value $1$.`
* **C4.** `$SG(\s_n) = 0$ if $n$ is even. Also, $SG(\s_n) = 1$ if $n = 1, 3, 5$ and $\mathcal \s_n = 2$ if $n$ is odd and 7 or larger.`
* **C5.** `$SG(\mathcal{\ps}_n) = 2$ if $n = 1, 3$ and $SG(\mathcal{\ps}_n) = 1$ if $n \geq 6$ is even and $SG(\mathcal{\ps}_5) = 5$ and $SG(\mathcal{\ps}_n) = 3$ otherwise.`
* **C6.** `\label{conj:hookswap} A partition $\lambda$ has Conway pair $(0, 1)$ if and only if $\lambda$ is an even-by-even hook.`
* **C7.** `A meld of losing thick hooks is losing.`
* **C8.** `\label{conj:bloby} Let $\lambda$ and $\Gamma_{r,c}$ be losing and suppose $\lambda^*=\Gamma_{r,c}\lhd\lambda$ is defined. Then $\lambda^*\in \p$.`
* … and 1 more.

### arXiv:2604.14639 — Unimodality and log-concavity of generalized Glasby-Paseman sequences
Seok Hyun Byun, Svetlana Poznanović (2026-04-16). 0 conjecture environment(s) in the TeX source.


### arXiv:2610.36132 — Hankel determinants for the $(n+2)$-shifted $q$-metallic numbers
Guo-Niu Han, Emmanuel Pedon (2026-10-26). 0 conjecture environment(s) in the TeX source.


### arXiv:2608.08703 — Hankel Transform and $(α,β)$ Somos-4 Sequences
Feihu Liu, Ying Wang, Zihao Zhang (2026-08-09). 0 conjecture environment(s) in the TeX source.


### arXiv:2602.24239 — Diamond Determinants and Somos Sequences
Nikolai Beluhov (2026-02-27). 6 conjecture environment(s) in the TeX source.

* **C1.** `\label{gre} Let $s$ be a Gale-Robinson sequence with nonzero terms, of a primitive type and even order $n$ with $n = 2m + 2$. Then the matrix $s \times s$ is of diamond rank at most $2^m$.`
* **C2.** `\label{gro} Let $s$ be a Gale-Robinson sequence with nonzero terms, of a primitive type and odd order $n$ with $n = 2m + 3$. Then the matrix $s \times s$ is of half-diamond rank at most $2^m$.`
* **C3.** `\label{cinv} For every proper type $\mathbf{n}$ of order $n$, it holds that $\dim \Omega_\boxtimes = \lfloor n/2 \rfloor$.`
* **C4.** `\label{cinve} In the setting of Conjecture \ref{cinv}: (a) Every element $\Phi$ of $\Omega_\boxtimes$ exhibits the symmetry $\Phi(x_0, x_1, \ldots, x_{n - 1}) = \Phi(x_{n - 1}, x_{n - 2}, \ldots, x_0)$; and (b) There exists a basis of $\Omega_\boxtimes$ which consists entirely of positive-coefficient polynomials when viewed over $\mathbb{Z}[\boldsymbol{\alpha}, \mathbf{x}]$.`
* **C5.** `\label{grte} Let $s$ and $t$ be two twinned Gale-Robinson sequences with nonzero terms, of a proper type and even order $n$ with $n = 2m + 2$. Then, generically, the matrix $s \times t$ is of diamond rank at most $2^m$.`
* **C6.** `\label{grto} Let $s$ and $t$ be two twinned Gale-Robinson sequences with nonzero terms, of a proper type and odd order $n$ with $n = 2m + 3$. Then, generically, the matrix $s \times t$ is of half-diamond rank at most $2^m$.`

### arXiv:2602.22129 — Enumeration of Nondegenerate $2 \times (k+1) \times k$ Hypermatrices
Brandon Koprowski, Joel Brewster Lewis (2026-02-25). 1 conjecture environment(s) in the TeX source.

* **C1.** `\label{mainConjecture} Fix a format $2\times (k+1) \times k$ and let $P = \brac{\lambda, \mu}$ be a plane partition such that $P \subseteq \Delta$. The number of nondegenerate $2\times (k+1) \times k$ hypermatrices over $\FF_q$ which respect $P$ is \[q^{k^2}(q-1)^{2k}[k+1 - \lambda_1]_q[k - \lambda_2]_q \cdots [2 - \lambda_{k}]_q \cdot [k - \mu_1]_q[(k-1) - \mu_2]_q \cdots[1 - \mu_k]_q.\]`

### arXiv:2502.06032 — A positivity conjecture for a quotient of $q$-binomial coefficients
Mona Gatzweiler, Christian Krattenthaler (2025-02-09). 0 conjecture environment(s) in the TeX source.


### arXiv:2502.21050 — Hankel determinants for convolution powers of Motzkin numbers
Ying Wang, Yingrui Zhang (2025-02-28). 3 conjecture environment(s) in the TeX source.

* **C1.** `For $ r \equiv 0 \pmod 3$, we have \begin{align*} H_{rn}(F(x,r))&=H_{rn+1}(F(x,r))=(-1)^{n}(n+1)^{r-1}.\\ H_{rn+2}(F(x,r))&+H_{rn-1}(F(x,r))= 3^2\alpha (n+1)^{r-1}. \end{align*} where $\lvert \alpha \rvert=\frac1{18}\,r \left( r-3 \right)$ .`
* **C2.** `For $ r \equiv 1 \ or \ 2 \pmod 3$, we have \begin{align*} H_{3rn}(F(x,r)) &= H_{3rn+1}(F(x,r)) = H_{3rn+r}(F(x,r)) = H_{3rn+r+1}(F(x,r))=\alpha.\\ H_{3rn+2r}(F(x,r))&=H_{3rn+2r+1}(F(x,r))= 0.\\ H_{3rn+2r-1}(F(x,r))&+H_{3rn+2r+2}(F(x,r))=\beta \left((2r)(n+1)\right)^{r-2}.\\ H_{3rn+2}(F(x,r))&+H_{3rn-1}(F(x,r))=\frac12\gamma r(r-3). \end{align*} where $\lvert \alpha \rvert=\lvert \beta \rvert=\lvert \gamma \rvert =1$.`
* **C3.** `For $ r \equiv 0 \pmod 3$, the shifted periodic of $H_n(F(x,r))$ is $r$. For $ r \equiv 1\ or\ 2 \pmod 3$, the shifted periodic of $H_n(F(x,r))$ is $3r$.`

### arXiv:2406.10987 — Inequalities for $k$-regular partitions
Bernhard Heim, Markus Neuhauser (2024-06-16). 1 conjecture environment(s) in the TeX source.

* **C1.** `Let $k \geq 2$. Let $N_k$ be the smallest integer, such that $\{p_k(n)\}_n$ is log-concave at every $n \geq N_k$. Then $N_k = N_{\infty}$ if $k\geq 30$.`

### arXiv:2406.13790 — The extended reverse ultra log-concavity of transposed Boros-Moll sequences
James J. Y. Zhao (2024-06-19). 6 conjecture environment(s) in the TeX source.

* **C1.** `\cite{Boros-Moll-2004}\label{conj:B-M-inflv} The Boros-Moll sequence $\{d_\ell(m)\}_{\ell=0}^m$ is $\infty$-log-concave.`
* **C2.** `The transposed Boros-Moll sequences $\{d_\ell(m)\}_{m\geq \ell}$ are $\infty$-strictly-log-concave for any $\ell\geq 3$.`
* **C3.** `For each $\ell\ge 1$, the sequence $\{d_\ell(m-1)d_\ell(m+1)/d_\ell^2(m)\}_{m\ge \ell+1}$ is log-concave.`
* **C4.** `For each $\ell\ge 0$, the sequence $\{d_\ell(m-1)d_\ell(m+1)/d_\ell^2(m)\}_{m\ge \ell+1}$ is extended reverse ultra log-concave.`
* **C5.** `\label{conj:hlc-hlv} For $m\ge 3$, the sequence $\{d_\ell(m)/d_{\ell-1}(m)\}_{1\le \ell \le \lfloor m/2 \rfloor+1}$ is log-convex, and the sequence $\{d_\ell(m)/d_{\ell-1}(m)\}_{\lfloor m/2 \rfloor\le \ell \le m}$ is log-concave. That is, \begin{align}\label{ineq:rlcv-l} \left(\frac{d_\ell(m)}{d_{\ell-1}(m)}\right)^2 < \left(\frac{d_{\ell-1}(m)}{d_{\ell-2}(m)}\right) \left(\frac{d_{\ell+1}(m)}{d_\ell(m)}\right) \end{align} holds for $2\le \ell \le \lfloor m/2 \rfloor$, and the inequalities in \eqref{ineq:rlcv-l} re …`
* **C6.** `The transposed Boros-Moll sequence $\{d_\ell(m)\}_{m\ge \ell}$ is infinitely log-monotonic for $\ell=0$. The sequence $\{d_\ell(m+1)/d_\ell(m)\}_{m\ge \ell}$ is infinite log-monotonic for each $\ell\ge 1$.`

### arXiv:2408.13346 — Partitions and elementary symmetric polynomials -- an experimental approach
Cristina Ballantine, George Beck, Mircea Merca (2024-08-23). 5 conjecture environment(s) in the TeX source.

* **C1.** `We have \begin{itemize} \item[(i)] $o_2(n)\geq d_2(n)$ for $n\geq 0$, \item[(ii)]$o_3(n)\geq d_3(n)$ for $n\geq 21$. \end{itemize}`
* **C2.** `\label{L_pre2} If $\lambda, \mu$ are partitions of $n$ for some $n\geq 0$, then $\pre_2(\lambda)=\pre_2(\mu)$ if and only if $\lambda=\mu$.`
* **C3.** `\label{L_prek} For $j\geq 3$, the mapping $\pre_j$ is injective on the set of partitions of $n$ of length greater than $j$.`
* **C4.** `\begin{align*}e_1p(n)^2 \geq e_1p(n-1)\cdot e_1p(n+1), & \text{ \ for } n>23; \\ e_2p(n)^2 \geq e_2p(n-1)\cdot e_2p(n+1),& \text{ \ for } n>17;\\ e_jp(n)^2 \geq e_jp(n-1)\cdot e_jp(n+1),& \text{ \ for } j>2, n>0. \end{align*}`
* **C5.** `For $j,n>0$, \begin{enumerate} \item[(i)] $\displaystyle{ \left( \frac{e_jp(n)}{p(n)} \right)^2 \geq \frac{e_jp(n-1)}{p(n-1)} \cdot \frac{e_jp(n+1)}{p(n+1)} }$,\\ \item[(ii)] $\displaystyle{ \left( \frac{e_jp(n)}{n\,p(n)} \right)^2 \geq \frac{e_jp(n-1)}{(n-1)\,p(n-1)} \cdot \frac{e_jp(n+1)}{(n+1)\,p(n+1)} }$. \end{enumerate}`

### arXiv:2604.17424 — Elementary symmetric polynomials and a potentially injective family of maps on partitions
Aman Devnani, Pramod Eyyunni (2026-04-19). 3 conjecture environment(s) in the TeX source.

* **C1.** `\label{conj12} For $n \geq 0$, the map $pre_2$ is injective on $\mathcal{P}(n)$.`
* **C2.** `\label{conj13} For $k \geq 3$, $pre_k$ is injective on the subset of $\mathcal{P}(n)$ consisting of partitions with at least $k$ parts.`
* **C3.** `\label{conj13new} For $k \geq 3$, $pre_k$ is injective on partitions of $n$ with $j$ parts for each $j \geq k+1$.`

### arXiv:2512.02324 — Parity results for the reciprocals of false theta functions
Jing Jin, Huan Xu, Olivia X. M. Yao (2025-12-02). 1 conjecture environment(s) in the TeX source.

* **C1.** `\cite{Keith} \label{conj-1} For $n\geq 0$, \begin{align} c_{13,1}(32n+23) &\;\equiv 0 \pmod 2,\label{1-8}\\[6pt] c_{13,1}(64n+63) &\;\equiv 0 \pmod 2,\label{1-9}\\[6pt] c_{13,1}(72n+j) &\;\equiv 0 \pmod 2,\label{1-10} \end{align} where $j\in\{15,21,39,69\}$.`

### arXiv:2512.04775 — Combinatorial proof of a result on generalized overcubic partitions and related conjectures
Suparno Ghoshal, Arijit Jana (2025-12-04). 1 conjecture environment(s) in the TeX source.

* **C1.** `\cite[Conjecture 7.3]{Das25}\label{p3conj1} For all $n\geq 0$ and $i\geq 1$, we have \begin{align*} \bar a_{3i+2}(9n+2)&\equiv 0 \pmod {6},\\ \bar a_{3i+2}(9n+5)&\equiv 0 \pmod{6},\\ \bar a_{3i+2}(9n+8)&\equiv 0 \pmod{6},\\ \bar a_{9i+5}(9n+3)&\equiv 0 \pmod{12},\\ \bar a_{9i+8}(9n+3)&\equiv 0 \pmod{12}.\\ \end{align*}`

### arXiv:2503.17571 — Finite analogs of partition bias related to hook length two and a variant of Sylvester's map
Alexander Berkovich, Aritram Dhar (2025-03-21). 3 conjecture environment(s) in the TeX source.

* **C1.** `(Ballantine et al. \cite[Conjecture $1.7$]{BBCFW23})\label{conj14} \quad\\ \begin{enumerate}[(i)] \item For every integer $t\ge 2$, there exists an integer $N_t$ such that for all $n > N_t$, we have $a_t(n)\ge b_t(n)$. Moreover, we conjecture the following values of $N_t$ for $2\le t\le 10$:\\ \begin{center} \begin{tabularx}{.8\textwidth} { | >{\centering\arraybackslash}X | >{\centering\arraybackslash}X | >{\centering\arraybackslash}X | >{\centering\arraybackslash}X | >{\centering\arraybackslash}X | >{\centering\ar …`
* **C2.** `\label{conj26} For every positive integer $L$, we have\\ \begin{align*} \sum\limits_{m\ge 0}\binom{m}{2}a_2(L,m,n)\ge \sum\limits_{m\ge 0}\binom{m}{2}b_2(L,m,n),\\ \end{align*} for all $n\ge 0$.\\`
* **C3.** `\label{conj41} For every positive integer $L$, we have\\ \begin{align*} \dfrac{1}{(q;q^2)_L}\left(\dfrac{1-q^{2L-2}}{1-q^2}+q^3\left[\begin{matrix}L-2\\2\end{matrix}\right]_{q^2}\right)\succeq q\left[\begin{matrix}L-2\\2\end{matrix}\right]_{q}(-q^3;q)_{L-4}.\numberthis\label{eq41}\\ \end{align*}`

### arXiv:2404.18517 — Distributions of statistics on separable permutations
Joanna N. Chen, Sergey Kitaev, Philip B. Zhang (2024-04-29). 3 conjecture environment(s) in the TeX source.

* **C1.** `\label{unimod-conj-rmax-sep} The distribution of separable permutations $\pi\in S_n$ with $\stat(\pi)=k$, $\stat\in\{\lmax,\rmax, \lmin, \rmin\}$, is unimodal with the peak being at $k=2$ for $n\geq 3$ (see Table~\ref{tab-rmax-sep}).`
* **C2.** `\label{conj-2} The distribution of $\rmax$ or $\lmin$ (resp., $\lmax$ or $\rmin$) on {\em irreducible} (resp., {\em reducible}) separable permutations is unimodal with the peak being at $k=3$ for $n\geq 5$ (see Table~\ref{tab-rmax-irr-sep}).`
* **C3.** `\label{conj-3} The distribution of $\rmax$ or $\lmin$ (resp., $\lmax$ or $\rmin$) on {\em reducible} (resp., {\em irreducible}) separable permutations is unimodal with the peak being at $k=1$ for $n\geq 1$ (see Table~\ref{tab-rmin-irr-sep}).`

### arXiv:2406.02971 — Maximal number of subword occurrences in a word
Wenjie Fang (2024-06-05). 6 conjecture environment(s) in the TeX source.

* **C1.** `\label{conj:limit-value-2} We have $L_2 > \logb(3/2)$.`
* **C2.** `For $k \geq 2$, let $w$ be a word of length $n \geq 1$ achieving the minimal subword entropy $\minswentropy^{(k)}(n)$. Then, except for a finite number of $n$, the longest run in $w$ has length $3$. Furthermore, the average run length converges when $n \to +\infty$.`
* **C3.** `There is a finite number of $n$ such that there is a binary word $w$ of length $n$ achieving the minimal subword entropy with several most frequent subwords.`
* **C4.** `\label{conj:periodic-subword-struct} For a given word $w$, there is a word $v$ such that, for all $m$ large enough, there is a most frequent subword of $w^m$ that takes the form $u \cdot v^r \cdot u'$, with $u$ and $u'$ of lengths bounded by $|v|$.`
* **C5.** `\label{conj:increasing-entropy} For fixed $k \geq 2$, there is a value $N$ such that the function $\minswentropy^{(k)}(n)/n$ is increasing for $n \geq N$.`
* **C6.** `\label{conj:half-length} For a given word $w$ of length at least $2$, there is at least one most frequent subword of $w$ with length at most $\lceil |w| / 2 \rceil$.`

### arXiv:2402.15886 — Extension of Bressoud's generalization of Borwein's conjecture and some exact results
Alexander Berkovich, Aritram Dhar (2024-02-24). 2 conjecture environment(s) in the TeX source.

* **C1.** `\label{conj13} Let $K$ be a positive integer and $N,M,\alpha K,\beta K$ be non-negative integers such that $1\le \alpha+\beta\le 2K-1$ (strict inequalities when $K = 2$) and $\beta-K\le N-M\le K-\alpha$. Then, $G(N,M;\alpha,\beta,K)$ is a polynomial in $q$ with non-negative coefficients.\\`
* **C2.** `\label{conj21} Let $K,i$ be positive integers such that $0 < i < K$ and $N,M,\alpha K,\beta K,\alpha i,\beta i$ be non-negative integers such that $1\le \alpha+\beta\le K-1$ (strict inequalities when $K = 4$ and $i=2$) and $\beta-i\le N-M\le K-\alpha-i$. Then, $D_{K,i}(N,M;\alpha,\beta)$, as defined in Theorem \ref{thm11}, is a polynomial in $q$ with non-negative coefficients.\\`

### arXiv:2401.09497 — Weakly Consecutive Sequences
Thomas Garrison, Chris Seiler, Andrew Knowles (2024-01-17). 6 conjecture environment(s) in the TeX source.

* **C1.** `\label{conj:main} Any WCS can be derived from one of the four starting sequences by applying the three rules. Moreover, for $k \geq 12$, the starting sequence and rule applications are unique up to reordering the rules.`
* **C2.** `\label{conj:1 at end} Any WCS either begins or ends with the term $1$ (in other words, $\sigma(1) = 1$ or $\sigma(k) = 1$).`
* **C3.** `\label{conj:m near end} Let $\sigma$ be a WCS. Then each $m < k/2$ appears as either one of the first $m$ or last $m$ terms of the sequence; in other words \[ \sigma^{-1}(m) \in [1,m] \cup [k-m+1,k]. \]`
* **C4.** `\label{conj:power of 2} The number $N(k)$ of WCS's is always a power of $2$ unless $k$ is of the form $2p - 1$ where $p$, $2p + 1$, and $p + 2$ are all primes. In the latter case, $N(k)$ is the sum of two powers of $2$.`
* **C5.** `If all WCS's are listed in order of increasing length $k$, then asymptotically 100\% of them are derived from Starting Sequence 4.`
* **C6.** `The value $N(k) = 2$ (minimal for $k > 1$) occurs for infinitely many $k$'s. Specifically, the density of $k$'s such that $N(k) = 2$ is given by \[ \delta = \prod_{p \text{ prime}} \left(1 - \frac{\log\left(1 + \frac{1}{p}\right)}{\log(p)}\right) \approx 0.19. \]`

### arXiv:2409.02789 — The entries of the Sinkhorn limit of an $m \times n$ matrix
Eric Rowland, Jason Wu (2024-09-04). 14 conjecture environment(s) in the TeX source.

* **C1.** `\label{general polynomial} Let $m \geq 1$ and $n \geq 1$. For every $m \times n$ matrix~$A$ with positive entries, the top-left entry $x$ of $\Sink(A)$ satisfies \[ \sum_{S \subseteq D(m, n)} \det\!\left(\tfrac{1}{m} \adj_S(m, n)\right) M(S) x^{\size{S}} = 0. \] In particular, $x$ is algebraic over the field generated by the entries of $A$, with degree at most $\binom{m + n - 2}{m - 1}$.`
* **C2.** `\label{polynomial form} Let $n \geq 1$. There exist integers $c_S(n)$, indexed by subsets $S \subseteq D(n, n)$, such that, for every positive $n \times n$ matrix~$A$, the top-left entry $x$ of $\Sink(A)$ satisfies \begin{equation}\label{sum over subsets} \sum_{k = 0}^{\binom{2 n - 2}{n - 1}} \left(\sum_{\substack{S \subseteq D(n, n) \\ \size{S} = k}} c_S(n) M(S)\right) x^k = 0 \end{equation} or, equivalently, \[ \sum_{S \subseteq D(n, n)} c_S(n) M(S) x^{\size{S}} = 0. \]`
* **C3.** `\label{symmetry} Let $n \geq 2$. For each $(R, C) \in D(n, n)$, define $ \kappa\!\left(\!\begin{smallmatrix} R\vphantom{\{} \\ C\vphantom{\{} \end{smallmatrix}\!\right) = \begin{smallmatrix} \{2, 3, \dots, n\} \setminus R \\ \{2, 3, \dots, n\} \setminus C \end{smallmatrix} $, and extend $\kappa$ to subsets of $D(n, n)$ by \[ \kappa\!\left(\!\begin{smallmatrix} R_1 & R_2 & \cdots & R_k\vphantom{\{} \\ C_1 & C_2 & \cdots & C_k\vphantom{\{} \end{smallmatrix}\!\right) = \left\{ \kappa\!\left(\!\begin{smallmatrix} R_1\v …`
* **C4.** `\label{binomial sum} For each $n \geq 1$, we have \[ \sum_{S \subseteq D(n, n)} c_S(n) x^{\size{S}} = (x - 1)^{\binom{2 n - 2}{n - 1}}. \]`
* **C5.** `\label{4 by 4 polynomial - large} There are $17920$ rational numbers \begin{align*} c_{D(4, 4)}(4) &= 1 \\ &\;\; \vdots \\ c_{\subsetmatrix{ \{\} \\ \{\} }}\!(4) &= -4 \\ c_{\subsetmatrix{ \{2\} \\ \{2\} }}\!(4) &= -2 \\ c_{\subsetmatrix{ \{2, 3\} \\ \{2, 3\} }}\!(4) &= 0 \\ c_{\subsetmatrix{ \{2, 3, 4\} \\ \{2, 3, 4\} }}\!(4) &= 2 \\ c_{\{\}}(4) &= 1 \end{align*} such that the top-left entry $x$ of the Sinkhorn limit of a $4 \times 4$ matrix satisfies \[ \sum_{k = 0}^{20} \left(\sum_S c_S(4) \Sigma(S)\right) x^k = …`
* **C6.** `\label{general polynomial form} Let $m \geq 1$ and $n \geq 1$. There exist rational numbers $c_S(m, n)$, indexed by subsets $S \subseteq D(m, n)$, such that, for every positive $m \times n$ matrix~$A$, the top-left entry $x$ of $\Sink(A)$ satisfies \[ \sum_{S \subseteq D(m, n)} c_S(m, n) M(S) x^{\size{S}} = 0. \]`
* **C7.** `\label{coefficient for size-1 subsets} Let $m \geq 1$ and $n \geq 1$, and let $ S = \subsetmatrix{ R \vphantom{\{} \\ C \vphantom{\{} } \subseteq D(m, n) $. Then $m c_S(m, n) = \size{R} (m + n) - m n$.`
* **C8.** `\label{coefficient for size-2 subsets} Let $m \geq 1$ and $n \geq 1$, and let $ S = \subsetmatrix{ R_1 & R_2\vphantom{\{} \\ C_1 & C_2\vphantom{\{} } \subseteq D(m, n) $. \begin{itemize} \item If $(R_1, C_1)$ and $(R_2, C_2)$ form a type-1 link with $\size{R_1} + 1 = \size{R_2}$, then \[ m^2 c_S(m, n) = \big(\size{R_1} (m + n) + n - m n\big) \big(\size{R_1} (m + n) + m - m n\big). \] \item If $(R_1, C_1)$ and $(R_2, C_2)$ form a type-2 link with $R_1 = R_2$, then \[ m^2 c_S(m, n) = \big(\size{R_1} (m + n) - m - m n …`
* … and 6 more.

### arXiv:2508.00946 — Partitioning set $[n] = \{1, \dots, n\}$ into subsets of size at most $m$ such that all sums are powers of $m$
Vladimir Gurvich, Mariya Naumova (2025-07-31). 1 conjecture environment(s) in the TeX source.

* **C1.** `\label{cnj-main} A $3$-good partition of $[n]$ exists for every positive integer $n$.`

### arXiv:2511.19744 — Toda primes
Stephen McKean (2025-11-24). 3 conjecture environment(s) in the TeX source.

* **C1.** `\label{conj:at least two} If $n$ is a positive integer, then $t(n)\geq 2$.`
* **C2.** `\label{conj:main} Let $n$ be an odd, square-free multiple of 3. Assume that there exists $p\in\{5,7,13\}$ such that $p\mid n$, and that $r\nmid n$ for $r\in\{5,7,13\}-\{p\}$. Finally, assume there exists $q\in T(3p)-\{5,7,13\}$ such that $q\nmid n$. Then $t(n)\geq 4$.`
* **C3.** `\label{conj:general bernoulli} Let $d$ be a Bernoulli denominator. \begin{enumerate}[(i)] \item\label{conj:toda for bernoulli} If $F(d)=4a$ for some integer $a$, then $T(am)=T(a)$ whenever $D_{4am}=D_{4a}$. \item\label{conj:bernoulli inequality} If $D_{4am}\neq D_{4a}=d$, then $t(am)\geq t(a)$. \end{enumerate}`

### arXiv:2604.17847 — Sophie Germain Primes and the Totient of Fibonacci Numbers
Aradhya Goel (2026-04-20). 1 conjecture environment(s) in the TeX source.

* **C1.** `\label{conj:conv} If $S(q) \neq \emptyset$ for an odd prime $q$, then there exists a prime $p \equiv 1 \pmod{q}$ with $z(p) \mid \pi(q)$.`

### arXiv:2609.01562 — Points of maximal traffic on a grid with obstruction
Juan Gil, Zhenni Liang, Ayodeji Odetola, Michael Weiner (2026-09-01). 1 conjecture environment(s) in the TeX source.

* **C1.** `\label{conj:496} For every $n\ge 496$ and every obstruction $B$ on the antidiagonal $x+y=n$, the maximum of $f_B$ is attained at $(1,1)$ and $(n-1,n-1)$. Equivalently, $\rho(n)<1$ for all $n\ge 496$.`

### arXiv:2609.11072 — Rational Approximations for Reciprocals of Multiple Zeta Values and Trivariate Cauchy Numbers
Ce Xu, Jianqiang Zhao (2026-09-10). 0 conjecture environment(s) in the TeX source.


### arXiv:2607.08366 — Minimum modulus for the unique multiset-sum problem
José A. R. Fonollosa (2026-07-09). 1 conjecture environment(s) in the TeX source.

* **C1.** `\label{conj:global} For every $n \ge 2$ and every $N < 2^{\,n} - 2^{\lfloor \log_2 n \rfloor}$, no set of $n$ residues is valid mod $N$; that is, the super-increasing set attains the least valid modulus over all size-$n$ sets.`

### arXiv:2412.10214 — Thron-type continued fractions (T-fractions) for some classes of increasing trees
Veronica Bitonti, Bishal Deb, Alan D. Sokal (2024-12-13). 1 conjecture environment(s) in the TeX source.

* **C1.** `[Trivariate symmetries on vincular patterns] \hfill\break \label{conj.trivariate} $\!\!\!$We have the relations \begin{subeqnarray} P_n (1, q, r, s) & = & P_n (1, q, s, r) \\[1mm] P_n (p, 1, r, s) & = & P_n (p, 1, s, r) \\[1mm] P_n (p, q, 1, s) & = & P_n (q, p, 1, s) \\[1mm] P_n (p, q, r, 1) & = & P_n (q, p, r, 1) \label{eq.conj.trivariate} \end{subeqnarray}`

### arXiv:2510.22385 — Two $t$-analogues of the tree inversion enumerator
Sam Hopkins (2025-10-25). 2 conjecture environment(s) in the TeX source.

* **C1.** `\label{conj:des} $\widetilde{I}_n(q,t) = \sum_{\pi \in \mathrm{PF}(n)}q^{\mathrm{cosum}(a)}t^{\mathrm{des}(\mathrm{oc}(\pi))}$.`
* **C2.** `For $n\geq 1$, write \[I_n(-1,t)=\sum_{k=0}^{\lfloor (n-1)/2\rfloor}B(n,k)t^{n-1-k} \textrm{ and } \widetilde{I}_n(-1,t)=\sum_{k=0}^{n-2}C(n,k)t^{k+1}.\] Then: \begin{enumerate} \item $B(n,k)$ is the number of ``\emph{simsun permutations}'' in $S_{n-1}$ with $k$ descents (see, e.g.,~\cite{chow2011counting} and \url{https://oeis.org/A094503}); \item $C(n,k)$ is the number of alternating permutations in $S_n$ with $k$ ``\emph{big returns}'' (see, e.g.,~\cite{petersen2025zigzag} and \url{https://oeis.org/A205497}). \e …`

### arXiv:2609.14295 — The Jensen--Pólya program and inequalities for arithmetic sequences
Koustav Banerjee, Kathrin Bringmann, Larry Rolen (2026-09-13). 14 conjecture environment(s) in the TeX source.

* **C1.** `\label{conj0} {\it For $k\ge 2$, $\overline{p}_k(n)$ satisfies \eqref{log-concave} and \eqref{HigherOrderTuranDef} for $n\gg 1$.} %$m\in\{2,3\}$, $J_{\overline{p}_k(n)}^{m,n}$ are hyperbolic for $n\!\gg\!1$.`
* **C2.** `\label{conj1} {\it For $k\in \mathbb{N}_{\ge 3}$ and $m\in \mathbb{N}$, the Jensen polynomials $J^{m,n}_{b_k}$ are hyperbolic for $n\!\gg\!1$.} %sufficiently large.`
* **C3.** `{\rm{(\cite[Conjecture 3.11]{LZ})}}\label{LZconj} {\it For $k\in \mathbb{N}_{\ge 2}$, $p_k(n)$ and $\overline{p}_k(n)$ satisfy \eqref{Briggsineq} for $n\!\gg\!1$.}`
* **C4.** `\label{Banconj} {\it We have},\footnote{The factor $n^{-\frac{3r}{2}}$ was missing in \cite{Ban}.} {\it as} $n\to \infty$, \[ L_r(p)(n)\sim \frac{(2r-1)!!\pi^r }{2^{\frac r2+5} 3^{\frac r2+1}}\frac{e^{2\pi\sqrt{\frac{2}3 n}}}{n^{\frac{3r}{2}+2}}. \]`
* **C5.** `\label{Mconj1} {\it We have, as} $n\to \infty$, \[ L_r\left(p_{24}\right)(n)\sim (2r-1)!!2^{r-2}\pi^r\frac{e^{8\pi\sqrt n}}{n^{\frac{3r+27}2}}. \]`
* **C6.** `\label{brokenLagconj1} {\it For $k\in\{1,2\}$ and $1\le r\le 14$, $\{b_k(n)\}_{n\ge N_{b_k}(r)}$ satisfies}\footnote{The values of $N_{b_k}(r)$ for $k\in\{1,2\}$, $1\le r\le 14$ were stated in \cite[Conjecture 6.2]{Yang1}.} \eqref{Lagineq}.`
* **C7.** `\label{Yangconj} {\it For $1\le r\le 5$, $\{\Delta^{[r]}(p)(n)\}_{n\ge T_p(r)}$ and $\{\Delta^{[r]}(\overline{p})(n)\}_{n\ge T_{\overline{p}}(r)}$ satisfy}\footnote{The values of $T_p(r)$ and $T_{\overline{p}}(r)$ were given in \cite[Table 5]{Yang}.} \eqref{HigherOrderTuranDef}.`
* **C8.** `\label{Yangconj1} {\it For $1\le m\le 11$ and $1\le r\le 5$, $\{\Delta^{[r]}(p)(n)\}_{n\ge L_p(r,m)}$ and}\ \ $\{\Delta^{[r]}(\overline{p})(n)\}_{n\ge L_{\overline{p}}(r,m)}$ {\it satisfy}\footnote{The values for $(L_{p}(r,m))_{\substack{1\le m\le 11\\1\le r\le 5}}$ and $(L_{\overline{p}}(r,m))_{\substack{1\le m\le 11\\1\le r\le 5}}$ were stated in \cite[Table 1 and 2]{Yang}.} \eqref{Lagineq}.`
* … and 6 more.

### arXiv:2609.05358 — Purely Periodic Three-move Subtraction Games
Hikaru Manabe (2026-09-04). 3 conjecture environment(s) in the TeX source.

* **C1.** `\label{conj:necessity} Let $a\ge2$, $\gcd(a,b,c)=1$ and $c\neq a+b$. Then $\G_S$ is purely periodic if and only if the angle is admissible for at least one candidate period, that is, if and only if \[ \rho\ \in\ \Th_{a+b}\cup\Th_{c+a}\cup\Th_{c+b} , \] and in that case the least period is the least admissible one, \[ p=\min\big\{P\in\{a+b,\,c+a,\,c+b\}:\ \rho\in\Th_P\big\}. \]`
* **C2.** `[period restriction]\label{conj:period} Let $S$ be non-additive with $a\ge2$ and $\gcd(a,b,c)=1$. If $\G_S$ is purely periodic then its least period satisfies $p\in\{a+b,\,c+a,\,c+b\}$.`
* **C3.** `[admissible angle]\label{conj:residue} Let $S$ be non-additive with $a\ge2$ and $\gcd(a,b,c)=1$, and suppose that $\G_S$ is purely periodic with least period $p$. If $p=c+a$ then $\rho\in\Th_{c+a}$, and if $p=c+b$ then $\rho\in\Th_{c+b}$.`

### arXiv:2609.16843 — On $t$-sumfree sequences
Daan van Berkel, Wieb Bosma (2026-09-15). 3 conjecture environment(s) in the TeX source.

* **C1.** `For every $g> f\geq 1$ the 2-sumfree sequence $S_{f, g}$ is ultimately periodic.`
* **C2.** `\label{conj:except} For every $d\geq 3$ the 3-sumfree sequence $S_{1, d-1, 2d-1}$ is ultimately periodic, and, with the exception of $d=3$ or $d=5$, the preperiod length is $d+3$ and the period length is $d$ modulo $8d-1$.`
* **C3.** `[folklore]\label{conj:folk} Every $3$-sumfree sequence $S_{f,g,h}$ is ultimately periodic.`

### arXiv:2404.18237 — Torus Queen Independence
Kada Williams (2024-04-28). 2 conjecture environment(s) in the TeX source.

* **C1.** `There is an absolute constant $c'>0$ such that for all $n$, either $Q_d(n)=n^{d-1}$ or $Q_d(n)\le n^{d-2}(n-c')$.`
* **C2.** `$Q_d(n)=n^{d-1}$ if and only if $(n,(2^d-1)!)=1$.`

### arXiv:2402.02199 — Neighborly boxes and bipartite coverings; constructions and conjectures
Jarosław Grytczuk, Andrzej P. Kisielewicz, Krzysztof Przesławski (2024-02-03). 5 conjecture environment(s) in the TeX source.

* **C1.** `\label{Conjecture Limit} For every fixed integer $k\geqslant 1$, there exists a real number $\gamma_k$ such that \begin{equation} \lim_{d\rightarrow \infty}\frac{n(k,d)}{d^k}=\gamma_k. \end{equation}`
* **C2.** `\label{Conjecture Pascal} For every $1\leqslant k\leqslant d$, \begin{equation} n(k,d)\leqslant n(k-1,d-1)+n(k,d-1). \end{equation}`
* **C3.** `\label{Conjecture n(2,d)=b(d)} For each $d\geqslant 2$, $n(2,d)=b(d)$.`
* **C4.** `\label{Conjecture bp} For every $n\geqslant 2$, we have $\bp_2(K_n)=c(n)$.`
* **C5.** `There exists a strategy in the splitting game resulting in a $k$-neighborly code of maximum size, that is, $\score(k,d)=n(k,d)$, for every $1\leqslant k\leqslant d$.`

### arXiv:2512.24221 — Grundy Values of Node-Kayles on Level-Regular Trees
Nuttanon Songsuwan (2025-12-30). 3 conjecture environment(s) in the TeX source.

* **C1.** `\label{conj:eventual} For all $a_1,\dots,a_{i-1}$ and $a_{i+1},\dots,a_d$, the sequence \[ n\mapsto\G\left(\T{a_1,\dots,a_{i-1},n,a_{i+1},\dots,a_d}\right) \] is eventually periodic with period $2$.`
* **C2.** `\label{conj:lep} Let $H$ be a level-regular rooted tree. Then for every graph $G$ and every vertex $u$ of $G$, the sequence $k\mapsto\G\left(G_u\langle H;k\rangle\right)$ is eventually periodic with period $2$.`
* **C3.** `\label{conj:runs} For every finite sequence $s$, the sequences $h\mapsto\G\left(\T{1^h,s}\right)$ and $h\mapsto\G\left(\T{s,1^h}\right)$ are eventually periodic, with period dividing $34$.`

### arXiv:2401.01385 — General Berndt-Type Integrals and Series Associated with Jacobi Elliptic Functions
Ce Xu, Jianqiang Zhao (2024-01-01). 1 conjecture environment(s) in the TeX source.

* **C1.** `Set $\Gamma=\Gamma(1/4)$. Then we have \begin{align*} \int_0^\infty \frac{x^9 dx}{(\cos x+\cosh x)^6} =&\,\frac{-63}{5\cdot 2^{10}} + \frac{1071\Gamma^{8}}{5^2\cdot 2^{13}\pi^{2}} -\frac{21\Gamma^{8}}{2^{12}\pi^{3}} +\frac{63\Gamma^{8}}{2^{16}\pi^{4}} -\frac{21\Gamma^{16}}{5^3\cdot 2^{13}\pi^{4}} + \frac{3\Gamma^{16}}{5\cdot 2^{13}\pi^{5}}\\ &\,-\frac{161\Gamma^{16}}{5\cdot 2^{19}\pi^{6}} +\frac{21\Gamma^{16}}{2^{21}\pi^{8}} +\frac{\Gamma^{24}}{5^2\cdot 2^{19}\pi^{8}} - \frac{\Gamma^{24}}{3\cdot 2^{20}\pi^{9}} + \f …`

### arXiv:2509.10023 — Some conjectures of Schlosser and Zhou on sign patterns of the coefficients of infinite products
Bing He, Linpei Li (2025-09-12). 6 conjecture environment(s) in the TeX source.

* **C1.** `\label{conj:S} Let $\delta$ be a real number satisfying \[ 0.227998127341\cdots\approx\frac{9-\sqrt{73}}{2}\leq\delta\leq1\quad\text{or}\quad2\leq\delta\leq3. \] Then the series $A^{\delta}(q),B^{\delta}(q),C^{\delta}(q)$ appearing in the dissection \[ G_{3}(q)^{\delta}=A^{\delta}(q^{3})-qB^{\delta}(q^{3})-q^{2}C^{\delta}(q^{3}) \] are power series in $q$ with non-negative real coefficients.`
* **C2.** `\label{conj:1} Let \[ Q_{5}(q):=\frac{(q,q^{4};q^{5})_{\infty}}{(q^{2},q^{3};q^{5})_{\infty}}. \] Then the $q$-series coefficients of $Q_{5}(q)^{\delta}$ exhibit the sign pattern $+-+--$ for \[ 1\leq\delta\leq\frac{\sqrt{97}-5}{2}\approx2.424428900898\cdots \] For \[ 2.571366313289\cdots\approx\alpha\leq\delta\leq4, \] they exhibit the sign pattern $+-+-+.$ Here $\alpha$ is the unique real root of the polynomial \[ x^{7}+35x^{6}+7x^{5}-6055x^{4}-14336x^{3}+104300x^{2}-184752x+282240 \] that satisfies $2<\alpha<3.$  …`
* **C3.** `\label{conj:2} Let \[ Q_{6}(q):=(q,q^{5};q^{6})_{\infty}. \] Then the $q$-series coefficients of $Q_{6}(q)^{\delta}$ exhibit the sign pattern $(+-)^{3}$ for all $\delta\geq3.$`
* **C4.** `\label{conj:3} Let \[ Q_{8}(q):=\frac{(q,q^{7};q^{8})_{\infty}}{(q^{3},q^{5};q^{8})_{\infty}}. \] Then the q-series coefficients of $Q_{8}(q)^{\delta}$ exhibit for $\delta=2$ the length $16$ sign pattern $+-++-+--+--+-++-.$ For \[ 2.664479110226972\cdots\approx\beta\leq\delta\leq4 \] they exhibit the sign pattern $+-++-+--.$ Here $\beta$ is the unique real root of the polynomial \[ \ensuremath{\begin{aligned} & x^{12}-90x^{11}+1457x^{10}+30486x^{9}-537081x^{8}\\ & +1892346x^{7}-3683653x^{6}-837509646x^{5}+774767020 …`
* **C5.** `\label{conj:4} Let \[ Q_{10}(q):=\frac{(q,q^{9};q^{10})_{\infty}}{(q^{3},q^{7};q^{10})_{\infty}}. \] Then the $q$-series coefficients of $Q_{10}(q)^{\delta}$ exhibit for $\delta=1$ the sign pattern $+-++--+--+.$ And for $\delta=-1$ the sign pattern $++++-----+.$`
* **C6.** `\label{conj:5} Let \[ Q_{12}(q):=\frac{(q,q^{11};q^{12})_{\infty}}{(q^{5},q^{7};q^{12})_{\infty}}. \] Then the $q$-series coefficients of $Q_{12}(q)^{\delta}$ exhibit for $\delta=1$ the sign pattern $+-+\,0-+-+-0+-,$ for $2\leq\delta\leq3$ they exhibit the sign pattern $+-+--+-+-++-.$ For $\delta=-1$ they exhibit the sign pattern $+++++\,0-----0,$ and for $-1<\delta<0$ the sign pattern $+++++------+.$`

### arXiv:2611.28840 — Congruences modulo powers of $2$ and $3$ for Lin's function $\text{PDO}_t(n)$
Manjil P. Saikia (2026-11-21). 1 conjecture environment(s) in the TeX source.

* **C1.** `\label{conj:sharp} For any $k\geq3$, there exists an odd integer $\lambda_k$ such that \eqref{eq:sharp} holds for any $n\geq0$.`

### arXiv:2609.11813 — Generalized Frobenius Partitions Modulo Powers of $2$
Manjil P. Saikia (2026-09-10). 1 conjecture environment(s) in the TeX source.

* **C1.** `[{\cite[Conjecture~5.1]{DasNathSarma2026}}] \label{conj:DNS} For every $n\geq0$, \begin{equation} \cphi_{18}(30n+19)\equiv0\pmod{16},\qquad \cphi_{18}(30n+25)\equiv0\pmod{16}. \label{eq:conjecture} \end{equation}`

### arXiv:2509.08801 — Arithmetic properties of partition functions introduced by Pushpa and Vasuki
Hemjyoti Nath, Manjil P. Saikia (2025-09-10). 1 conjecture environment(s) in the TeX source.

* **C1.** `\cite[Conjecture 7.1]{DasappaChannabasavayyaKeerthana}\label{conjd} For all $n\geq 0$ and $\alpha\geq 1$, we have \[ K(7^\alpha n+7^\alpha -2)\equiv 0 \pmod{7^\alpha}. \]`

### arXiv:2508.00074 — On simultaneous $(s, s+t, s+2t, \dots)$-core partitions
William Keith, Rishi Nath, James Sellers (2025-07-31). 2 conjecture environment(s) in the TeX source.

* **C1.** `[Fayers \cite{Fayers}]\label{fayersconj} Suppose $t \geq 1$. Then there is a monic polynomial $f_t(s)$ of degree $t-1$ with non-negative integer coefficients such that for any $s \geq 1$ coprime to $t$, the number of $(s,s+t,s+2t,\dots)$-cores is $$\frac{2^{s-t}f_t(s)}{t!}.$$ The constant term of $f_t(s)$ is $(2^t-1)(t-1)!$.`
* **C2.** `[Fayers, \cite{Fayers}]\label{fayersconj2} Consider the same environment as Conjecture \ref{fayersconj}, and suppose $t \geq 2$. Then $f_t(s)$ is divisible by $s+t+(-1)^t$.`

### arXiv:2404.15716 — Parity of the coefficients of certain eta-quotients, III: two special classes
William J. Keith, Fabrizio Zanello (2024-04-24). 5 conjecture environment(s) in the TeX source.

* **C1.** `\label{nontrivfivehalves} If $2^\alpha < 5t_0/3$, $\alpha > 0$, then both $$\frac{f_{5t}^3}{f_1} {\ }{\ }\textnormal{and}{\ }{\ } \frac{f_{t/2} f_{5t/2}}{f_1}$$ have density $1/2$.`
* **C2.** `[\cite{KZ}, Conjecture 4]\label{mainconj} Let $F(q)=\sum_{n\ge 0} c(n)q^n$ be an eta-quotient, shifted by a suitable power of $q$ so powers are integral, and denote by $\delta_F$ the odd density of its coefficients $c(n)$. We have: i) For any $F$, $\delta_F$ exists and satisfies $\delta_F\le 1/2$. ii) If $\delta_F= 1/2$, then for any nonnegative integer-valued polynomial $P$ of positive degree, the odd density of $c(P(m))$ is $1/2$. (In particular, $c(Am+B)$ has odd density $1/2$ for all arithmetic progressions $Am …`
* **C3.** `Let $t= t_0 2^{\alpha}$ be even, with $t_0$ odd, $t_0 > 3 \cdot 2^{\alpha}$ (i.e., $G_t(q)$ does not satisfy the conditions of CMSZ). Then $\delta^{(t)} = 1/2$.`
* **C4.** `Let $t \equiv 1 \pmod{4}$, $t \geq 29$. Then there exists no arithmetic progression $An+B$, $A \neq 0$, such that $g_t(An+B)$ is always even. Moreover (which is equivalent under Conjecture \ref{mainconj}), we have $\delta^{(t)} = 1/2$.`
* **C5.** `\label{evenmodulus} Let $t$ be odd. If $g_t(An+B) \equiv 0 \pmod{2}$ for all $n$, then $A$ is even.`
