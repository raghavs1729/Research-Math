# Literature / novelty check, 27 Nov 2026 (krtw4)

Scope: KRTW (Kanade–Russell–Tsuchioka–Warnaar, arXiv:2404.03851v2) identities (ii), (iii), (iv) after
Conjecture 4.9 for k >= 3; Conjecture 2.2 (general z); Conjecture 4.9; Primc's D_3^{(2)} conjecture, level >= 3.
Builds on `../rr-identities/PRIOR-ART.md`, section "KRTW second round" (26 Nov 2026), which is not redone here.

Legend: [read TeX] = arXiv source downloaded and read around the relevant places; [abstract only] = abs page /
API metadata only; [not accessible] = refused or unavailable.

Status: DONE (27 Nov 2026).

**Bottom line.** No proof, and no claimed proof, of KRTW (ii), (iii) or (iv) for any k ≥ 3, of Conj. 2.2 for general z
(k ≥ 2), of Conj. 4.9 (k ≥ 3), or of Primc's $D_3^{(2)}$ / three-row conjecture at level ≥ 3 was found in: arXiv
(everything in math.CO/NT/QA/RT submitted 19–25 Nov 2026 per the API and announced through Thu 26 Nov per the listing
pages, plus targeted queries), OpenAlex and Semantic Scholar citation lists of KRTW (8 citing works each; all but one
already in the prior record), GitHub code/issue/PR search, and the two trackers. The only
new KRTW-citing items found are Šikić 2404.01791 ($C_2^{(1)}$, 2024) and Dombos–Jang 2611.17443 (bib entry only, not
actually cited); neither is relevant. KRTW's own text (v2 of 20 Apr 2026; the journal version is Selecta 2026) says:
the four identities are "conjectural" (derived from Conj. 4.9 + CMPP/Primc conjectures), k = 2 of Conj. 4.9 is Thm 4.10,
Conj. 2.2 is stated as known for k = 0, 1 only, the $n=2$ (three-row) case of Primc's conjecture "appears to still be
open", and the $w,z$-multisums for the other weights ($1\le a\le k-2$) are an "open problem". Verbatim statements:
§1 below (1b = the four identities).

## 1. Verbatim statements from arXiv:2404.03851v2 [read TeX]

Source: `https://arxiv.org/e-print/2404.03851v2` (single file `KRTW.tex`, 145 511 bytes; downloaded to the scratchpad,
not committed). Abstract page: v1 5 Apr 2024, v2 20 Apr 2026; journal ref "Selecta Mathematica, New Series 32 (2026),
no. 5, Paper No. 100"; 45 pages. Numbering below is by the LaTeX counters of the arXiv TeX (theorem counter shared by
all environments, reset per section; the journal numbering was not checked):
Conj. 2.2 = `Con_Shun`; Conj. 4.1 = `Con_Cn1` (CMPP); Conj. 4.2 = `Con_Dn2` (Primc); Conj. 4.7 = `Con_Cn1Dn2-qseries`;
Prop. 4.8 (unlabelled); Conj. 4.9 = `Con_Shun2`; Thm. 4.10 = `Thm_kis2`. Macros: `\la`=\lambda, `\La`=\Lambda,
`\Lap`=\Uplambda (upright Lambda, used for the fundamental weights), `\C`=\mathscr{C}, `\abs{x}`=|x|.

### 1a. Conjecture 2.2 (`Con_Shun`), TeX lines 1021–1043, with the sentence before it and the remark after it

```latex
From \eqref{Eq_PS-2lambda} (which is a special case of \eqref{Eq_P-spec}),
\[
\sum_{\substack{\la \\[1pt] \la_1\leq k}}
(zq)^{\abs{\la}} P_{2\la}(1,q,q^2,\dots;q)=
\sum_{r_1,\dots,r_k\geq 0}\,
\prod_{i=1}^k \frac{z^{r_i} q^{r_i^2}}{(q;q)_{r_i-r_{i+1}}}.
\]
Conjecturally, this generalises to Hall--Littlewood functions of base $q^2$.

\begin{conjecture}\label{Con_Shun}
For $k$ a nonnegative integer,
\begin{equation}\label{Eq_q2}
\sum_{\substack{\la \\[1pt] \la_1\leq k}}
(zq)^{\abs{\la}} P_{2\la}(1,q,q^2,\dots;q^2)
=\sum_{\substack{r_1,\dots,r_k\geq 0 \\[1pt] s_1,\dots,s_k\geq 0}}
\, \prod_{i=1}^k \frac{z^{r_i+s_i} q^{(r_i+s_i)^2+s_i^2+s_i}}
{(q;q)_{r_i-r_{i+1}}(q^2;q^2)_{s_i-s_{i-1}}},
\end{equation}
where $r_{k+1}=s_0:=0$.
\end{conjecture}

Apart from the trivial case $k=0$ this also holds for $k=1$ by
\eqref{Eq_GOW27} for $n=1$, $\delta=0$ and $(r,r_1)\mapsto (r_1+s_1,s_1)$.
```

(`Eq_GOW27` = Proposition 2.1 (`Prop_GOW`): $P_{(2^r)}(1,q,q^2,\dots;q^{2n+\delta})=\sum_{r_1,\dots,r_n\ge0}
q^{r^2-r+r_1^2+\dots+r_n^2+r_1+\dots+r_n}/((q;q)_{r-r_1}(q;q)_{r_1-r_2}\cdots(q;q)_{r_{n-1}-r_n}(q^{2-\delta};q^{2-\delta})_{r_n})$,
"stated without proof in [GOW16, Equation (2.7)]", proved in KRTW by the Bailey lemma.)

### 1b. Conjecture 4.9 (`Con_Shun2`) and the four "conjectural Andrews–Gordon-type identities", TeX lines 2435–2501

```latex
By \eqref{Eq_q2}, the $n=2$ instance of \eqref{Eq_D} admits
an alternative expression as a multisum.
We conjecture that similar such multisums hold for
$\mathcal{D}^{(2)}_{k\Lap_1}(z,q)$ and
$\mathcal{D}^{(2)}_{\La_0+(k-1)\Lap_1}(z,q)$.

\begin{conjecture}\label{Con_Shun2}
For $k$ a nonnegative integer,
\begin{align*}
\mathcal{D}^{(2)}_{k\Lap_0}(z,q)&=
\sum_{\substack{r_1,\dots,r_k\geq 0 \\[1pt] s_1,\dots,s_k\geq 0}}\,
\prod_{i=1}^k \frac{z^{r_i+s_i} q^{(r_i+s_i)^2+s_i^2+r_i+2s_i}}
{(q;q)_{r_i-r_{i+1}}(q^2;q^2)_{s_i-s_{i-1}}}, \\
\mathcal{D}^{(2)}_{k\Lap_1}(z,q)&=
\sum_{\substack{r_1,\dots,r_k\geq 0 \\[1pt] s_1,\dots,s_k\geq 0}}\,
\prod_{i=1}^k \frac{z^{r_i+s_i} q^{(r_i+s_i)^2+s_i^2}}
{(q;q)_{r_i-r_{i+1}}(q^2;q^2)_{s_i-s_{i-1}}}, \\
\mathcal{D}^{(2)}_{\Lap_0+(k-1)\Lap_1}(z,q)&=
\sum_{\substack{r_1,\dots,r_k\geq 0 \\[1pt] s_1,\dots,s_k\geq 0}}\,
\Omega_{r_1,\dots,r_k}^{s_1,\dots,s_k}(q)
\prod_{i=1}^k \frac{z^{r_i+s_i} q^{(r_i+s_i)^2+s_i^2}}
{(q;q)_{r_i-r_{i+1}}(q^2;q^2)_{s_i-s_{i-1}}},
\end{align*}
where
\[
\Omega_{r_1,\dots,r_k}^{s_1,\dots,s_k}(q):=
\sum_{i=1}^{k-1}q^{r_i+2s_i}(1-q^{2s_{i+1}-2s_i})+q^{r_k+2s_k}
\]
and $r_{k+1}=s_0:=0$.
\end{conjecture}

Combining this with \eqref{Eq_CDn2} and Conjectures \ref{Con_Cn1}
and \ref{Con_Dn2} gives the following conjectural Andrews--Gordon-type 
identities:
\[
\sum_{\substack{r_1,\dots,r_k\geq 0 \\[1pt] s_1,\dots,s_k\geq 0}}\,
\prod_{i=1}^k \frac{q^{(r_i+s_i)^2+s_i^2+s_i}}
{(q;q)_{r_i-r_{i+1}}(q^2;q^2)_{s_i-s_{i-1}}}
=\frac{(q,q^{k+1},q^{k+2};q^{k+2})_{\infty}}
{(q;q^2)_{\infty}(q;q)_{\infty}},
\]
\begin{align*}
\sum_{\substack{r_1,\dots,r_k\geq 0 \\[1pt] s_1,\dots,s_k\geq 0}} &\,
\prod_{i=1}^k \frac{q^{(r_i+s_i)^2+s_i^2+r_i+2s_i}}
{(q;q)_{r_i-r_{i+1}}(q^2;q^2)_{s_i-s_{i-1}}} \\
&\qquad=
\frac{(q,q^2,q^{2k+2},q^{2k+3},q^{2k+4},q^{2k+4};q^{2k+4})_{\infty}}
{(q^2;q^2)_{\infty}(q;q)_{\infty}},
\end{align*}
\begin{align*}
\sum_{\substack{r_1,\dots,r_k\geq 0 \\[1pt] s_1,\dots,s_k\geq 0}} &\,
\prod_{i=1}^k \frac{q^{(r_i+s_i)^2+s_i^2}}
{(q;q)_{r_i-r_{i+1}}(q^2;q^2)_{s_i-s_{i-1}}} \\
&\qquad=
\frac{(q^{k+1},q^{k+2},q^{k+2},q^{k+3},q^{2k+4},q^{2k+4};q^{2k+4})_{\infty}}
{(q^2;q^2)_{\infty}(q;q)_{\infty}}
\end{align*}
and
\begin{align*}
\sum_{\substack{r_1,\dots,r_k\geq 0 \\[1pt] s_1,\dots,s_k\geq 0}}\, 
\Omega_{r_1,\dots,r_k}^{s_1,\dots,s_k}(q) &
\prod_{i=1}^k \frac{q^{(r_i+s_i)^2+s_i^2}}
{(q;q)_{r_i-r_{i+1}}(q^2;q^2)_{s_i-s_{i-1}}} \\
&\qquad=
\frac{(q^k,q^{k+1},q^{k+3},q^{k+4},q^{2k+4},q^{2k+4};q^{2k+4})_{\infty}}
{(q^2;q^2)_{\infty}(q;q)_{\infty}}.
\end{align*}
```

Checks against the brief: (i)–(iv) of the task statement match the print exactly, including the denominator of (iv),
which is printed as $(q^2;q^2)_\infty(q;q)_\infty$ (the brief wrote "(...)"), and $\Omega$, whose $i$-th term is
$q^{r_i+2s_i}(1-q^{2s_{i+1}-2s_i})$ for $1\le i\le k-1$ plus $q^{r_k+2s_k}$. Note the first display ((i)) uses
$q^{\dots+s_i}$, not $q^{\dots+r_i+2s_i}$. KRTW derive all four from Conj. 4.9 together with (Eq_CDn2) and Conjs 4.1, 4.2;
independently of Primc, (i) ⇔ Conj. 2.2 at $z=1$ and (ii) ⇔ Conj. 2.2 at $z=q$ (the HL sides are GOW Thm 1.2, $n=1$, and
Thm 1.3, $n=2$; see the 26 Nov record), while (iii), (iv) ⇔ lines 2, 3 of Conj. 4.9 at $z=1$ **plus** Primc's Conj. 4.2
($n=2$) for the weights $k\Uplambda_1$, $\Uplambda_0+(k-1)\Uplambda_1$ — there is no HL form for those weights
(Section 5 of KRTW, quoted below).
In the intro (TeX lines 633–653) the authors single out (i) and (iii): (i) "is a $q$-series identity for the principally
specialised character of the $\mathrm{A}^{(1)}_1$-standard module $L(k\La_0)$ very different from those in
\cite{Andrews10,CSXY22,Feigin09,KY13,Stembridge90,Tsuchioka23}", (iii) "is a $q$-series identity related to a non-standard
specialisation of the $\mathrm{D}^{(2)}_3$-standard module $L(k\La_1)$". "The close resemblance of these formulas to the
Andrews--Gordon identities \eqref{Eq_AG} seems quite remarkable."

### 1c. Remarks on which cases are known (verbatim)

After the display (TeX lines 2503–2535):
```latex
Conjecture~\ref{Con_Shun2} for $k=2$ can be completed to a full
set of weights, which is provable using the functional equations
for $\mathcal{D}^{(2)}_{\La}$.

\begin{theorem}\label{Thm_kis2}
We have
\begin{align*}
\mathcal{D}^{(2)}_{2\Lap_0}(z,q)&=
\sum_{\substack{r_1,r_2\geq 0 \\[1pt] s_1,s_2\geq 0}}\,
\prod_{i=1}^2 \frac{z^{r_i+s_i}q^{(r_i+s_i)^2+s_i^2+r_i+2s_i}}
{(q;q)_{r_i-r_{i+1}}(q^2;q^2)_{s_i-s_{i-1}}}, \\
\mathcal{D}^{(2)}_{2\Lap_1}(z,q)&=
\sum_{\substack{r_1,r_2\geq 0 \\[1pt] s_1,s_2\geq 0}}\,
\prod_{i=1}^2 \frac{z^{r_i+s_i} q^{(r_i+s_i)^2+s_i^2}}
{(q;q)_{r_i-r_{i+1}}(q^2;q^2)_{s_i-s_{i-1}}}, \\
\mathcal{D}^{(2)}_{\Lap_0+\Lap_1}(z,q)&=
\sum_{\substack{r_1,r_2\geq 0 \\[1pt] s_1,s_2\geq 0}}\,
\prod_{i=1}^2 \frac{z^{r_i+s_i}
q^{(r_i+s_i)^2+s_i^2+\delta_{i,2}(r_i+2s_i)}}
{(q;q)_{r_i-r_{i+1}}(q^2;q^2)_{s_i-s_{i-1}}}\,
\Big(1+zq^{2+\sum_i(r_i+2s_i)}\Big), \\
\mathcal{D}^{(2)}_{\Lap_0+\Lap_2}(z,q)&=
\sum_{\substack{r_1,r_2\geq 0 \\[1pt] s_1,s_2\geq 0}}\,
\prod_{i=1}^2 \frac{z^{r_i+s_i}
q^{(r_i+s_i)^2+s_i^2+r_i+2\delta_{i,2}s_i}}
{(q;q)_{r_i-r_{i+1}}(q^2;q^2)_{s_i-s_{i-1}}}\,
\Big(1+zq^{2+\sum_i(r_i+2s_i)}\Big),
\end{align*}
where $r_3=s_0:=0$.
\end{theorem}

The expression for $\mathcal{D}^{(2)}_{\Lap_0+\Lap_1}$ agrees with
the one given in Conjecture~\ref{Con_Shun2}.
```
(The proof, TeX 2556–2748, solves a w-deformed system of four functional equations by "the computer-assisted procedure
of \cite{KR23}".)

Conjecture 4.7 (`Con_Cn1Dn2-qseries`) and Proposition 4.8, TeX lines 2364–2433 (statement and the claims; proof abridged):
```latex
\begin{conjecture}\label{Con_Cn1Dn2-qseries}
For $k$ a nonnegative integer,
\begin{subequations}
\begin{align}\label{Eq_C}
\mathcal{C}^{(n)}_{k\Lap_0}(z,q)=
\sum_{\la\in\C_{k,0^n}} z^{l(\la)} q^{\abs{\la}}
&=\sum_{\substack{\la \\[1pt] \la_1\leq k}}
(zq)^{\abs{\la}} P_{2\la}\big(1,q,q^2,\dots;q^{2n}\big) 
\intertext{and}
\mathcal{D}^{(n)}_{k\Lap_0}(z,q)=
\sum_{\la\in\D_{k,0^n}} z^{l(\la)} q^{\abs{\la}}
&=\sum_{\substack{\la \\[1pt] \la_1\leq k}}
(zq^2)^{\abs{\la}} P_{2\la}\big(1,q,q^2,\dots;q^{2n-2}\big),
\label{Eq_D}
\end{align}
where in \eqref{Eq_C} it is assumed that $n\geq 0$ and
in \eqref{Eq_D} that $n\geq 1$.
\end{subequations}
\end{conjecture}

\begin{proposition}
Equation \eqref{Eq_C} holds for $n=0$, \eqref{Eq_D} holds
for $n=1$, and both equations hold for $k=1$.
Moreover, \eqref{Eq_C} holds for $z=1$.
\end{proposition}
[...]
Finally, the $z=1$ case of \eqref{Eq_C} follows by combining
\cite[Theorem 1.2]{GOW16} with the recent
Primc--Trup\v{c}evi\'c proof of Conjecture~\ref{Con_Cn1}
for $k_0=k$ and $k_2=\dots=k_n=0$.
```

Primc's conjecture (Conj. 4.2, `Con_Dn2`, TeX 1888–1920; product as given by KRTW) and its status:
```latex
\begin{conjecture}[Primc {\cite[Conjecture 2.1]{Primc24}}]\label{Con_Dn2}
For $k$ a nonnegative integer, $n$ a positive integer and
$k_0,\dots,k_n$ nonnegative integers such that $k_0+\dots+k_n=k$,
\begin{align*}
&\sum_{\la\in\D_{k_0,\dots,k_n}} q^{\abs{\la}} \\
&\qquad= \frac{(q^{2k+2n};q^{2k+2n})_{\infty}^n}{(q^2;q^2)_{\infty}
(q;q)_{\infty}^{n-1}}
\prod_{1\leq i<j\leq n}
\theta\big(q^{\la_i-\la_j-i+j},q^{\la_i+\la_j+2n-i-j+1};q^{2k+2n}\big),
\end{align*}
where $\la_i:=k_i+\dots+k_n$ for $1\leq i\leq n$.
\end{conjecture}

Primc does not actually state what the product on the right should be.
[...]
Unlike Conjecture~\ref{Con_Cn1} however, the case of three-row
frequency arrays (i.e., $n=2$) appears to still be open.
```

End of Section 4 (TeX 2772–2823), on w-deformations and the missing members of the family:
```latex
We do not yet know what the appropriate set of $w$-deformed equations
is for $k\geq 3$ with the exception of the two equations
\[
\mathcal{D}^{(2)}_{k\Lap_0}(z,w)=\sum_{i=0}^k 
(zq^2)^i\,\mathcal{D}^{(2)}_{i\Lap_0+(k-i)\Lap_1}(zq^2,wq^2)
\]
and
\[
\mathcal{D}^{(2)}_{k\Lap_1}(z,w)-
\mathcal{D}^{(2)}_{\Lap_0+(k-1)\Lap_1}(z,w)
=(zq)^k \sum_{i=0}^k \Big(\frac{wq}{z}\Big)^i\,
\mathcal{D}^{(2)}_{i\Lap_0+(k-i)\Lap_1}(zq^2,wq^2).
\]
[...]
A reasonable guess is thus that
\begin{align*}
\mathcal{D}^{(2)}_{k\Lap_0}(z/q,w/q^2,q)&=
\mathcal{D}^{(2)}_{k\Lap_1}(z,w,q) \\
&=\sum_{\substack{r_1,\dots,r_k\geq 0 \\[1pt] s_1,\dots,s_k\geq 0}}\,
\prod_{i=1}^k \frac{z^{r_i} w^{s_i} q^{(r_i+s_i)^2+s_i^2}}
{(q;q)_{r_i-r_{i+1}}(q^2;q^2)_{s_i-s_{i-1}}}, \\
\mathcal{D}^{(2)}_{\La_0+(k-1)\Lap_1}(z,w,q)
&=\sum_{\substack{r_1,\dots,r_k\geq 0 \\[1pt] s_1,\dots,s_k\geq 0}}
\Omega_{r_1,\dots,r_k}^{s_1,\dots,s_k}(q)
\prod_{i=1}^k \frac{z^{r_i} w^{s_i} q^{(r_i+s_i)^2+s_i^2}}
{(q;q)_{r_i-r_{i+1}}(q^2;q^2)_{s_i-s_{i-1}}},
\end{align*}
where $r_{k+1}=s_0:=0$.
[...]
The above series thus have the structure of two interwoven copies
of the Andrews--Gordon multisums.
It is an open problem to find the $w,z$-generalisations
of the Andrews--Gordon multisums \eqref{Eq_AG} for $1\leq a\leq k-2$.
Of course, the correct such sum should equate to
\[
\frac{(q^{a+1},q^{a+2},q^{2k-a+2},q^{2k-a+3},
q^{2k+4},q^{2k+4};q^{2k+4})_{\infty}}
{(q^2;q^2)_{\infty}(q;q)_{\infty}}
\]
when $z=w=1$.
```
(Matching products with the general-$a$ product: (ii) is the case $a=0$, (iii) the case $a=k$, (iv) the case $a=k-1$ of
a conjectural $(k+1)$-member family (indexing as in Eq_AG2, TeX 537–541, where $a$ = number of $r_i$ *without* a linear
term: $q^{r_1^2+\dots+r_k^2+r_{a+1}+\dots+r_k}$, product $(q^{a+1},q^{2k-a+2},q^{2k+3};q^{2k+3})_\infty/(q)_\infty$); members
$1\le a\le k-2$ have no conjectured multisum ("open problem"). In this "reasonable guess" the weight $k\Uplambda_0$ is
obtained from $k\Uplambda_1$ by $(z,w)\mapsto(zq,wq^2)$, i.e. (ii) = (iii) with the extra factor $q^{\sum(r_i+2s_i)}$.)

Start of Section 5 (TeX 2839–2844): "For the weights $\La=k\Lap_0$ and $\La=k\Lap_n$ explicit formulas for these generating
functions in terms of Hall--Littlewood symmetric functions are proposed in Conjectures~\eqref{Con_A2n2-qseries} and
\eqref{Con_Cn1Dn2-qseries}. No such expressions in terms of Hall--Littlewood functions seem to exist for other weights."

Not stated anywhere in the TeX: no numerical-verification range is given for Conj. 2.2 / 4.9 / (i)–(iv) (grep for
"verified", "checked", "numerically", "computer", "tested": only the computer-assisted *proof* of Thm 4.10).

## 2. New arXiv papers, 2026-11-20 .. 2026-11-27, and citation lists

**Systematic window sweep [abstract only].** arXiv API, `cat:math.X AND submittedDate:[202611190000 TO 202611272359]` for
X = CO (318 entries), NT (193), QA (37), RT (76) (567 unique; latest submission date indexed: 25 Nov), plus the listing
pages `/list/math.X/new` (announcement of Thu 26 Nov 2026) and `/list/math.X/pastweek` (20, 23, 24, 25, 26 Nov). Titles and
abstracts filtered for Rogers/Ramanujan/partition/Bailey/Hall–Littlewood/q-series/Primc/Capparelli/Andrews/Gordon/
Bressoud/multisum/Nahm/Kanade/Warnaar/Tsuchioka/Meurman/fermionic/Macdonald/q-difference/hypergeometric/... Every
q-series-type hit (abstract page opened via the API):

| arXiv | date | content | relevance |
|---|---|---|---|
| 2611.32425 Capparelli–Meurman–Primc, "Three combinatorial identities of the Rogers-Ramanujan type" | 24 Nov | three identities motivated by $A_2^{(2)}$, proved by the Takigiku–Tsuchioka automata method | already in the prior record; $A_2^{(2)}$, not $D_3^{(2)}$ |
| 2611.26719 Stanton, "Partial sums of Rogers-Ramanujan functions" | 19 Nov | positive multisums on base $q^5$ for the RR sum and its partial sums; a quintic basic hypergeometric identity | none; TeX grepped: no Kanade/Warnaar/Primc/Capparelli/2404.03851 [read TeX, grep only] |
| 2611.25684 Hughes, "A sharp bounded bijection for tight cylindric partitions" | 19 Nov | tight cylindric partitions, Warnaar's $Q_{n,c}$ | none; cites Kanade–Russell only, not KRTW [read TeX, grep only] |
| 2611.29412 Li–Milas, "Principal subspaces of affine vertex algebras at admissible levels" | 22 Nov | classification, product character formulas ($\mathfrak{sp}_4$ at $-1/2$, $G_2$ at $-5/3$, $C_3$) | none; no Kanade/Primc/Capparelli/2404.03851 in TeX [read TeX, grep only] |
| 2611.29740 Hakimoglu-Brown, "Two-term Rogers dilogarithm identities over cubic fields" | 22 Nov | dilogarithm identities | none [abstract only] |

Other targeted API queries (sorted by submission date, newest first): `all:Primc` (47 hits; newest 2611.32425),
`all:Capparelli` (55; newest 2611.32425), `abs:"Andrews-Gordon"` (69; newest in the window: none — newest is
2609.20567), `au:Kanade` (newest math paper 2606.19234), `au:Warnaar` (newest 2607.08607, sandpile groups). Nothing
posted 20–27 Nov 2026 beyond 2611.32425.

**OpenAlex [metadata].** KRTW = W4394591349 (journal version: Selecta Math. (N.S.) 32 (2026), no. 5, Paper 100,
DOI 10.1007/s00029-026-01198-y), `cited_by_count` = 8: Russell 2306.16251 (SIGMA 22 (2026) 046), Warnaar 2511.17034
(SIGMA 22 (2026) 062), Lau–Ono 2608.05480, Kurşungöz 2605.13096, Šikić 2606.19994, Bhattacharya–Dai–Hou–Tong
2610.41155, Bartholdi 2611.04020, Ding 2611.10654 — all in the prior record (2611.32425 not yet indexed by OpenAlex).
Works citing Russell 2306.16251 (10) add only Kanade–Russell "Tight cylindric partitions" (Proc. AMS, DOI
10.1090/proc/17685), Butorac–Kožić–Meurman–Primc (J. Algebra 2024), Layne–Marshall–Sadowski–Shambaugh (Ramanujan J.
2024), Šikić (Eur. J. Math. 2024), Kanade 2606.19234 — none on $D_3^{(2)}$ or the two-chain sums.

**Semantic Scholar** `paper/arXiv:2404.03851/citations`: HTTP 429 on the first attempt, 200 on the second: 8 citing
papers — Ding 2611.10654, Bartholdi 2611.04020, Bhattacharya–Dai–Hou–Tong 2610.41155, Lau–Ono 2608.05480, Šikić
2606.19994, Warnaar 2511.17034, Russell 2306.16251, and **Šikić 2404.01791** ("Combinatorial relations among relations of
$C_2^{(1)}$-standard modules for higher levels", Eur. J. Math. 2024; not in the prior record, but it is about
$C_2^{(1)}$ relations among relations, not $D_3^{(2)}$ or multisums [title/metadata only]).

**Found via GitHub (see §3): arXiv:2611.17443**, Dombos–Jang, "$A_1^{(1)}$-Grounded partitions at levels 1 and 2, Part II:
affine crystal structure" (12 Nov 2026, math.CO/math.RT) [abstract + TeX grep]. Its `references.bib` contains a KRTW
entry, but the TeX has no `\cite{KRTW}` (and no `\nocite`), so KRTW is not actually cited. Content: $A_1^{(1)}$ levels 1
and 2 only (crystal structure on grounded partitions; new $q$-series identities at levels 1, 2). Not relevant to
(ii)–(iv). (Part I is 2508.02664v3, 29 Jun 2026.)

Also noticed (A_2^{(1)}, not $D_3^{(2)}$): arXiv:2610.40079, Chern–Nguyen–Tran–Wakhare, "The bivariate generating
function for cylindric partitions of rank three" (29 Oct 2026) [abstract only] — finishes conjectures of Warnaar and
Kanade–Russell on rank-3 cylindric partitions.

**Conclusion of §2: no paper posted 20–27 Nov 2026 (or found in any citation list) proves or even addresses (ii), (iii),
(iv) for k ≥ 3, Conj. 2.2 for general z, Conj. 4.9, or Primc's $D_3^{(2)}$ conjecture at level ≥ 3.**

## 3. GitHub searches and trackers

`gh api` search endpoints only (27 Nov 2026):

| query | result |
|---|---|
| `search/code q='2404.03851'` and `q='"2404.03851"'` | HTTP 408 "Request timed out" both times (not retried further) |
| `search/issues q='2404.03851 is:issue'` / `is:pull-request` | 0 / 0 |
| `search/issues q='Kanade Warnaar conjecture is:issue'` / `is:pull-request` | 0 / 0 |
| `search/issues q='Capparelli Meurman Primc is:issue'` / `is:pull-request` | 0 / 0 |
| `search/code q='Capparelli Meurman Primc'` | 46 files: arXiv RSS/abstract mirrors (ehijano/rss_fetch, megaelius/CAI, IDBIDO/caim, Epc333/Uni, Xsardas1000/Search, LeeMath-lzu/daily-arXiv-ai-enhanced-QA), lenis2000/homepage `_arxiv/2024-04-05-2404.03851.md` (abstract card), MatthewCRussell homepage, dujella.github.io, jikokennobun/My-Reserch-Project (literature inbox, OpenAlex card for KRTW, no assessment), AlexenderSokolov/cross-theory-formal-database (see below) |
| `search/code q='Kanade Warnaar conjecture'` | 47 files: same mirrors plus rrsdb/rrsdb.github.io (Rogers–Ramanujan identity database: only bibliography entries GOW16, KLRS15, KR15 — no KRTW entry), tomdif/RamanujanLean `docs/CYLINDRIC_RESEARCH.md` (rank ≥ 4 cylindric-partition computations, refers to 2610.40079; "No novelty claim"), qkwang98/Arxivs (arXiv source mirrors: `2611.17443_source/references.bib` — led to Dombos–Jang, §2; `2611.18198_source` cites Kanade–Linshaw and Warnaar's $_1\psi_1$ note only), chineseschrodingercat/elvia (`studies/math2025/evidence/math2025-0254`: an audit of Warnaar's $A_2$ Bailey tree paper arXiv:2303.09069, JEMS — reports a "Major" gap in "the printed proof of the central Bailey kernel identity" with a repair; unrelated to KRTW, but noted since that paper is a technique source) |

AlexenderSokolov/cross-theory-formal-database: `editable-corpus/sources/189546/source0.tex` is a mirror of the TeX of
Warnaar 2511.17034 (SIGMA style, cites "KRTW24 ... Selecta Math. (N.S.), to appear"); its final section only says that
"at least conjecturally, the GOW identities of Corollaries 1–3 admit interpretations as identities for the generating
function of restricted sets of coloured partitions" and asks for interpretations of the other 14 — nothing on (ii)–(iv).
`items/1259_*` / `sources/1259/*` mirror a paper whose "second and third authors" are Kanade and Russell (intro surveys
$A_9^{(2)}$, $A_2^{(2)}$, $A_2^{(1)}$ identities; only the matched fragments were seen, paper not identified).

Trackers (raw.githubusercontent.com, HEAD): `1337isnot1337/solved-conjectures` README.md (94 table rows) and the full
catalogue CONTENTS.md (415 lines) — no entry for Primc, Capparelli, Kanade, Warnaar, Tsuchioka, Andrews–Gordon,
Rogers–Ramanujan, Bressoud, Bailey, multisum, q-series (the only Hall–Littlewood hit is an unrelated Chen–Khare–Sahi
Macdonald entry). `openai/math` README.md and CONTENTS.md (649 KB) — no such entry (the only "Ramanujan" hits concern
the generalized Ramanujan conjecture for automorphic forms).

## 4. Techniques (Kanade–Lovejoy, Dousse–Jang, Warnaar, Chern–Li–Stanton–Xue–Yee)

Summary: **none of the four papers contains a multisum with two coupled chains** (a base-$q$ chain $r_1\ge\dots\ge r_k$
and a base-$q^2$ chain $s_1\le\dots\le s_k$ coupled through $q^{(r_i+s_i)^2}$), and none has a rank-two / two-variable
Bailey step. A grep of all four TeX files for $(r_i+s_i)^2$-type exponents and for lines containing both $(q;q)_{\cdot}$
and $(q^2;q^2)_{\cdot}$ factors found nothing except one Warnaar example (below), which has a different shape.

* **Chern–Li–Stanton–Xue–Yee**, "The Ariki–Koike algebras and Rogers–Ramanujan type partitions", arXiv:2209.07713,
  J. Algebraic Combin. 60 (2024), no. 2, 491–540 [read TeX, §§1, 5; theorem numbers not checked, TeX labels given].
  Theorems `th:gen-sum`, `thm2.1`: for
  $m\ge1$, $0\le a\le m$, the generating function of Kleshchev multipartitions $\Lambda^{a,m}$ ($q=Q_1=\dots=Q_a=-1$,
  $Q_{a+1}=\dots=Q_m=1$) is
  $\sum_{N_1,\dots,N_m\ge0} q^{\sum_i\binom{N_i+1}{2}}/(q;q)_{N_1}\prod_{i}\left[{N_i+\delta_{m-a,i}\atop N_{i+1}}\right]
  =(q^{a+1},q^{m+1-a},q^{m+2};q^{m+2})_\infty/((q;q)_\infty(q;q^2)_\infty)$ — one base-$q$ chain of strict-partition
  weights linked by $q$-binomials. **Its product side at $a=0$ (or $a=m$) with $m=k$ is exactly the product side of KRTW
  (i)** (KRTW list [CSXY22] among the "very different" sum sides for $L(k\Lambda_0)$ of $\mathrm{A}^{(1)}_1$). Technique:
  a symmetry $a\leftrightarrow m-a$ (`th:sym-a`), a $q$-binomial multisum transformation removing the "+1"
  (`th:2a-general`), then
  Andrews' parity identity (2010, eq. (1.9)) and Kim–Yee (2013). Also (§4) bivariate (2-residue) generating functions
  for $m=2$ via Berkovich–Uncu and Zeng, giving base-$q^2$ double sums such as
  $\sum_{r,s} q^{r^2+s^2+r+s}(q^2;q^2)_{r+s+1}/((q^2;q^2)_r^2(q^2;q^2)_s(q^2;q^2)_{s+1})=(-q^2;q^2)_\infty/(q^2;q^2)_\infty$.
* **Kanade–Lovejoy**, "Bailey pairs and an identity of Chern–Li–Stanton–Xue–Yee", arXiv:2410.15408v2, SIGMA 21 (2025),
  021 [read TeX, §§1–3, 6; TeX labels given]. Theorem `thm:main`: for any Bailey pair relative to $q$ with $\alpha_n=\frac{1-q^{2n+1}}{1-q}f_n$,
  $q^nf_n=-q^{-n-1}f_{-n-1}$: $\sum q^{\sum\binom{n_i+1}{2}}(q^2;q^2)_{n_1+\delta_{a,0}}\beta_{n_1+\delta_{a,0}}/(q)_{n_m}
  \prod[{n_{i+1}+\delta_{a,i}\atop n_i}]=\frac{(-q)_\infty}{(q)_\infty}\sum_n q^{m\binom n2+an}f_n$; the CLSXY identity is
  the seed $\beta_n=1/(q^2;q^2)_n$ (Slater E(3)). Further families with $(q^2;q^2)$-denominators and $q^2$-binomials
  (`q^2thm`: $\sum q^{\sum(n_i^2+n_i)}/((q^2;q^2)_{n_m}(-q;q^2)_{n_1+\delta_{a,0}})\prod[\cdot]_{q^2}=
  (-q^2;q^2)_\infty(q^{2a+1},q^{2m-2a+2},q^{2m+3};q^{2m+3})_\infty/(q^2;q^2)_\infty$), false-theta companions, and dilated
  $\mathscr S_{m,a}$ (differences are products, `thm:dilationex1`). Technique: classical Bailey lemma,
  Bailey lattice (Agarwal–Andrews–Bressoud), and Lovejoy's lemma (2022) turning a pair relative to $x$ into one relative
  to $xq$ with $\beta'_n=(1-q^{n+1})\beta_{n+1}$ — this is what produces the shifted $q$-binomials
  $[{n_{i+1}+1\atop n_i}]$. Single chain; the base is mixed only in the seed. "w-deformed $\mathscr D$" in its concluding
  remarks refers to Kanade–Milas–Russell (Nahm sums with double poles, arXiv:2109.01909), not to KRTW's $\mathcal D^{(2)}$.
* **Dousse–Jang**, "Andrews–Gordon type identities with parity restrictions through particle motion",
  arXiv:2603.05300v1 (5 Mar 2026) [read TeX, §§1, 6; TeX labels given]. Theorems `thm:main`, `thm:even1`, `thm:even2`:
  single **base-$q^2$** chains
  $\sum_{s_1\ge\dots\ge s_k}q^{\sum s_i^2+\text{linear}}/((q^2;q^2)_{s_1-s_2}\cdots(q^2;q^2)_{s_k})$ = (finite sums of)
  products $(-q^2;q^2)_\infty(q^{2c},q^{2k+4-2c},q^{2k+4};q^{2k+4})_\infty/(q^2;q^2)_\infty$ or
  $(-q;q^2)_\infty(\dots;q^{2k+4})/(q^2;q^2)_\infty$ (Stanton-type generalisations of Andrews' and Kim–Yee's parity
  identities), with partition interpretations (odd/even parts appear an even number of times). Corollary `cor:2` gives
  a new proof of CLSXY by differencing $\mathrm{AK}_{a,k}-\mathrm{AK}_{a+1,k}$ and $q\to q^{1/2}$. Technique: Warnaar's
  particle-motion bijection (as developed by Dousse–Jang–Jouhet–Konan). No two-chain sums; no mention of KRTW/Primc.
  Open problems: combinatorial proofs of two sum-sum identities, binomial (Stanton-type) extensions.
* **Warnaar**, "Affine Jacobi–Trudi identities and $q,t$-Rogers–Ramanujan identities", arXiv:2511.17034v2, SIGMA 22
  (2026), 062 [read TeX: introduction and "Open problems"; prior pass also read it]. Conjectural affine/Hall–Littlewood dual Jacobi–Trudi
  formulas for $\mathrm{B}_n$, $\mathrm{C}_n$, $\mathrm{BC}_n$ HL polynomials of rectangular shape $(k^n)$ (`Con_JTB`,
  `Con_JTC`, `Con_JTBC`; special cases proved), implying **proved** $q,t$-identities (`Thm_qtRR`–`Thm_qtRR4`) of the form
  $\sum_{\lambda\text{ even},\lambda_1\le2k}q^{(\sigma+1)|\lambda|/2}P_\lambda(1,q,q^2,\dots;t)=\frac{(p;p)^k_\infty}{(q;q)^k_\infty}
  \prod\theta(\dots;p)$, $p=tq^{2k+2}$ ($\mathrm{C}^{(1)}_k$; analogues for $\mathrm{A}^{(2)}_{2k}$, $\mathrm{D}^{(2)}_{k+1}$,
  $\mathrm{A}^{(2)}_{2k-1}$), proved independently of the conjectures (section labelled `Sec_Ismail`). Specialising $t=q^N$ gives
  the GOW identities via a level-rank duality; with $\lambda=2\mu$, $t=q^2$, $\sigma=0$ (resp. 1) the left side is KRTW's
  $\sum_{\mu_1\le k}q^{|\mu|}P_{2\mu}(1,q,\dots;q^2)$ (resp. $q^{2|\mu|}$), i.e. the HL sides of KRTW (i), (ii) at $z=1$.
  **No multisum (fermionic) side for $t=q^2$**; only $t=q$ gives Andrews–Gordon sums (via
  $P_{(2^r)}(1,q,\dots;q)=q^{r^2-r}/(q)_r$). Nearest shape found: section "Open problems", for $t=q$,
  $\sum_{r,s\ge0}q^{(r+s)^2+2s^2+\sigma(r+s)}/((q;q)_r(q;q)_{2s})=\frac{(q^{14};q^{14})_\infty}{(q;q)_\infty}
  \theta(q^{2-\sigma};q^{14})\theta(q^{10+2\sigma};q^{28})$ (Virasoro $(3,7)$; from $\sum_{r,s}q^{(\sigma+1)(r+s)}t^{s^2}
  P_{(2^r,1^{2s})}$) — two coupled variables, but $(q;q)_{2s}$ rather than $(q^2;q^2)_s$, and a single level. KRTW is
  cited only for its Appendix A.2 and in the open problems ("it would be extremely interesting to find a partition
  theoretic interpretation ... for the remaining 14 corollaries"). Technique: Jacobi–Trudi determinants, Macdonald /
  Lassalle–Schlosser formulas, bounded Littlewood identities, Watson quintuple product.

Additional technique leads met on the way (abstract or partial reading only):
* Bhattacharya–Dai–Hou–Tong, "Changing parameters in Hall–Littlewood functions", arXiv:2610.41155 (30 Oct 2026)
  [abstract + TeX grep]: transition coefficients $b_{\lambda\mu}(q,t)$ of $P_\lambda(t)$ in the basis $P_\mu(q)$ —
  conceptually the $t=q\to t=q^2$ change of base behind Conj. 2.2 — but only hook and two-column $\lambda$ are treated
  (for $2\lambda$ with $\lambda_1\le k$ this covers only $k=1$, already known); KRTW is used only for the expansion of
  $P_{(2^r,1^s)}(t)$ in $e_{r-i}e_{r+s+i}$ (KRTW p. 11).
* Bressoud–Ismail–Stanton, "Change of base in Bailey pairs", math/9909053 [title/authors via arXiv API only]; used with
  the standard Bailey lemma by Dousse–Lovejoy arXiv:2509.21230 ("we require some changes-of-base established by
  Bressoud–Ismail–Stanton") [abstract only]. This is the classical tool for chains that mix $(q;q)$ and $(q^2;q^2)$
  factors; neither paper treats a two-chain sum of KRTW type.
* Chen–Liu–Sun arXiv:2605.06248 "Double-sum Rogers–Ramanujan type identities" ($q$-Hermite / Chebyshev Bailey pairs)
  [abstract only] — not KRTW shape.

## 5. Hall–Littlewood at alphabets {±a q^i}

Brief pass; the prior record already covered the citation lists of Bartlett–Warnaar 1304.1602, GOW and Rains–Warnaar.
arXiv API (abstracts): `abs:"Hall-Littlewood" AND abs:"Rogers-Selberg"` → only KRTW; `abs:"Hall-Littlewood" AND
abs:"Rogers-Ramanujan"` → 12 papers (Warnaar 2511.17034, KRTW, Larson 1506.08313 and 1506.05063, Frechette–Locus
1407.2880, GOW 1401.7718, Bartlett–Warnaar 1304.1602, Lascoux 0709.2088, Warnaar math/0410592, Ishikawa–Jouhet–Zeng
math/0401252, Jouhet–Zeng math/0110165, Fulman math/9712236) — none evaluates at an alphabet $\{\pm aq^i\}$ or at
$(y,-y)$ [titles/abstracts only, except as noted]; `abs:"Hall-Littlewood" AND (roots of unity OR plethysm)` → 8 papers,
about the *parameter* $t$ at roots of unity, not the alphabet. Items looked at:

* Frechette–Locus, arXiv:1407.2880, Ann. Comb. 20 (2016) 345–360 [abstract + theorem list of TeX]: the GOW series
  $R_\nu(a,b;q)=\sum_{\lambda_1\le a}q^{c|\lambda|}P_{2\lambda}(1,q,q^2,\dots;q^{2b+d})$ (for $b=1$, $d=0$ this is the HL
  side of KRTW Conj. 2.2 at $z=1$); only coefficient recursions derived from the product side and congruences — no
  multisum.
* Lascoux, arXiv:0709.2088, "Adding ±1 to the argument of an Hall–Littlewood polynomial" [abstract only]: the shift
  $p_i\to p_i\pm1$ (adjoining / removing the letter 1), not a $\pm$ pair of letters; recovers an identity of Warnaar.
* Descouens, arXiv:math/0701221 [abstract only]: HL functions of rectangular shape at primitive roots of unity
  (parameter) as plethysms, via ribbon rigged configurations.
* Warnaar 2511.17034 [TeX grep]: the $\mathrm{BC}_n$ HL polynomials carry factors $(1-s_1x_i)(1-s_2x_i)/(1-x_i^2)$, so a
  $\pm$ pair can occur as the *parameters* $(s_1,s_2)=(a,-a)$; no alphabet $(y,-y)$. KRTW TeX: no $\pm$ alphabets.

Elementary remark (not from a source): since $p_r(X\cup -X)=(1+(-1)^r)p_r(X)$, any symmetric function satisfies
$f(X\cup-X)=(V_2f)(X^2)$ with the Verschiebung $V_2\colon p_{2r}\mapsto 2p_r,\ p_{2r+1}\mapsto0$, so the
$(y,-y)$-evaluation of the $C_2$ Rogers–Selberg sum is a one-letter evaluation of $V_2$ of it; no source computing $V_2$
of (modified) HL functions in a form usable here was found in this pass.

**Conclusion of §5: nothing found** (consistent with the 26 Nov record): no fermionic formula for
$\sum_{\lambda_1\le k}q^{|\lambda|}P'_{2\lambda}(y,-y;q)$, for $P'_{2\lambda}(x,-x;q)$, or for HL sums at $\{\pm aq^i\}$.
