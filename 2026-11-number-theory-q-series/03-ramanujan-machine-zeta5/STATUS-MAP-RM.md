# STATUS MAP — Ramanujan Machine conjectures (polynomial continued fractions)

Survey date: 2026-11-24. Scope: which publicly posted Ramanujan Machine (RM) conjectures are still
unproved, and which are proved (with citation). No proofs were attempted here.

Evidence labels: **[F]** = I read the source in full or in the relevant part; **[A]** = abstract /
snippet / second-hand citation only; **[I]** = my inference. "num-check" = I evaluated the continued
fraction numerically (mpmath, 40–60 digits) only to make sure the formula is transcribed correctly;
a num-check is never evidence of proof status.

Convention used throughout (the RM site's and Yamamoto's): `PCF[a_n, b_n] = a_0 + b_1/(a_1 + b_2/(a_2 + b_3/(a_3 + ...)))`,
`a_n` for n >= 0, `b_n` for n >= 1. Cohen's database uses `(a(n), b(n)) = a(0) + b(0)/(a(1) + b(1)/(a(2)+...))`,
so (after his initial terms) RM `a_n` = Cohen `a(n+1)`, RM `b_n` = Cohen `b(n)`.

---

## 0. Headline

1. **Every itemized conjecture on the four site lists (zeta(3), Catalan, pi^2, mixed-zeta) and the log 2
   list has a proof in the literature.** The site's "new and unproven" labels are stale: the PDFs date
   from 2020–2022 and the Results page text was never updated [F]. Counts in section 1.
2. Caveat on "proved": the proofs of the 22 Catalan items, 7 of the pi^2 items and the log 2 item rest on
   **one unrefereed arXiv preprint** (Yamamoto, arXiv:2403.09729v1, 2024; no journal reference found [A]).
   The RM group itself, in arXiv:2610.27799v1 (20 Oct 2026), still writes "Ramanujan Machine Catalan limits
   remain conjectural" (Table 3 caption) and cites the catalogue as "22 entries marked 'new and unproven'"
   (ref. [88]) [F]. That paper cites neither Yamamoto nor Shvets [F: full-text grep, zero hits]. So there is a disagreement
   between the proof literature and the group's own bookkeeping; I found no published objection to
   Yamamoto's proof, and the group's statement reads as not-updated rather than as a rebuttal [I].
3. The second zeta(3) conjecture (Domb, 12/(7 zeta(3))) is proved only in a single-author April 2026
   preprint (Shvets, arXiv:2604.06239v1). The value itself was tabulated as a "known" Apéry limit since
   2008 and the RM group's own Ramanujan Library paper (2024) lists it as "Known" [F].
4. What is genuinely open is listed in section 2: three unidentified zeta(5)-level PCFs from the
   PNAS/CMF paper, two e-conjectures from the Nature paper for which I located no proof, and the general
   (structural) conjectures of the group's papers.

---

## 1. Tables per source list

### 1.1 zeta(3) list — `https://ramanujanmachine.com/wp-content/uploads/2020/06/zeta3.pdf` [F, text + rendered image]

Site text: "two of which are still unproven". 5 entries: 3 "known", 2 "new and unproven".

| id | formula (num-check done for all) | site label | status today |
|---|---|---|---|
| zeta3 #1 | 1/zeta(3) = PCF[n^3+(n+1)^3, -n^6] | known | classical (Euler) |
| zeta3 #2 | 5/(2 zeta(3)) = PCF[2+n(n+2)(3n+4), 4n^6-2n^5] | known | classical (series sum (-1)^(n-1)/(n^3 C(2n,n)); Cohen DB label 5.5.18 [F]) |
| zeta3 #3 | 6/zeta(3) = PCF[(2n+1)(17n(n+1)+5), -n^6] | known | Apéry |
| zeta3 #4 ("Z1") | **8/(7 zeta(3)) = PCF[(2n+1)(3n(n+1)+1), -n^6]** | new and unproven | **PROVED / pre-existing.** (a) The CF itself is in Batut–Olivier, Sém. Th. Nombres Bordeaux 1979–80, exp. 23, §3.2.4 (7 zeta(3)/8 = 1/(p(0) - 1^6/(p(1) - 2^6/...)), p(n)=6n^3+9n^2+5n+1) [A, via Shvets]. (b) Equivalent Apéry limit lim B_n/s_n = 7 zeta(3)/32 for s_n = sum_k C(n,k)^2 C(2k,n)^2 (recurrence (n+1)^3 u_{n+1} = (2n+1)(12n^2+12n+4) u_n - 16 n^3 u_{n-1}): Golyshev, "Deresonating a Tate period", arXiv:0908.1458, Thm 2.4 (variety V16) + Prop. 3.3 + Cor. 3.4 [F]; second proof Golyshev–Kerr–Sasaki, J. London Math. Soc. 109 (2024) e12825, §5.5 [A, via Shvets]. (c) Explicit rederivation and the continuant bridge to the RM normalisation: Shvets, arXiv:2604.14219v2 (24 Apr 2026), Thm 1.1 and intro [F intro]. (d) Cohen DB label 5.5.20, reached by Apéry acceleration from 5.5.8 and 5.5.10 [F]; Cohen arXiv:2304.11727 §15 [F]. (e) RM group's Ramanujan Library paper arXiv:2412.12361v2, App. F table: C[-n^6/(36n^6-21n^4+7n^2-1)] = 8/(7 zeta(3)) marked "Known" [F]. |
| zeta3 #5 ("Z2") | **12/(7 zeta(3)) = PCF[(2n+1)(5n(n+1)+2), -16n^6]** | new and unproven | **PROVED (unrefereed preprint).** Shvets, arXiv:2604.06239v1 (5 Apr 2026): Thm 1.1 (Domb Apéry limit lim B_n/D_n = 7 zeta(3)/24, D_n = sum_k C(n,k)^2 C(2k,k) C(2n-2k,n-k)), Thm 1.3 (the RM CF), Prop. 2.3 (P_n = (n+1)!^3 D_{n+1}/2^(n+1), Q_n = (n+1)!^3 B_{n+1}/2^n) [F §§1–2; modular part §§3–6 not checked]. Earlier: value listed as a known Apéry limit in Almkvist–van Straten–Zudilin 2008 without proof [A, per Shvets]; used "experimentally" in Cohen arXiv:2304.11727 §13 and Cohen DB label 5.5.17 [F]; marked "Known" in arXiv:2412.12361v2 App. F (C[-16n^6/(100n^6-45n^4+21n^2-4)] = 6/(7 zeta(3))) [F]. I found no proof earlier than Shvets (checked Yang 2008 [F grep], Golyshev 2009 [F], Chamberland–Straub 2021 [F grep]). |

CLAIMED on the forum: "yaron" (Chun Wing Lai), posts of June 2026 "Algebraic Proof ... 8/(7zeta(3))" and
"... 12/(7zeta(3))" — I read both images [F]: they only rewrite the CF as its Euler series
`a_0 - sum (prod b_i)/(Q_n Q_{n-1})`; the value is assumed, not derived. **Not proofs.**

### 1.2 Catalan list — `.../2020/06/catalan.pdf` [F]

23 entries: #1 "known", #2–#23 "new and unproven" (22). All 22 are **PROVED** in S. Yamamoto,
"Proof and generalization of conjectures of Ramanujan Machine", arXiv:2403.09729v1 (13 Mar 2024), as
Corollaries 2.1–2.22 (statements pp. 12–14, proofs on the pages that follow) [F statements; F for the proofs
of Thm 2.1 and Thm 2.2; the corollary proofs only skimmed].

Parametric family (Yamamoto Thm 2.1, valid for Re(2γ-α-β)<1, extended by analytic continuation in
Lemma 2.3/Thm 2.2):
`PCF[Δ(n(n+α)(n+β)), -2n(n+α)(n+β)(n+γ)] = 4 / sum_{n>=0} (1/2)_n (γ+1)_n / [((α+1)/2)_{n+1} ((β+1)/2)_{n+1}]`,
where Δf(n)=f(n+1)-f(n), i.e. a_n = 3n^2 + (2α+2β+3)n + (α+1)(β+1). Thm 2.2 gives two contiguity relations
in α and γ for H(α,β;γ) := 1/PCF. (num-check of Thm 2.1 at three generic parameter triples: agreement to
the accuracy of the slowly converging series.) The same families appear, *without proof and only up to a
Möbius transformation of G*, in Cohen's database arXiv:2409.06086v5 at label 5.4.29 [F]:
family A: `a(n)=3n^2+(4v+4w-3)n+(2v-1)(2w-1), b(n)=-2n(n+u)(n+2v)(n+2w)` (= Yamamoto with α=2v, β=2w, γ=u);
family B: `a(n)=3n^2+(2u+2v+2w-3)n+uv+uw+vw-u-v-w+1, b(n)=-2n(n+u)(n+v)(n+w)`.

| id | value = PCF[a_n, b_n] (all num-checked) | Yamamoto | family parameters [I] |
|---|---|---|---|
| cat #1 | 6/(8G - pi*acosh 2) = PCF[10n^2+7n+2, -(2n-1)^4-(2n-1)^3] | (site: known) | — |
| cat #2 | 1/(2G) = PCF[3n^2+3n+1, -2n^4] | Cor 2.1 | (α,β,γ)=(0,0,0); Cohen 5.4.29 itself |
| cat #3 | 2/(2G-1) = PCF[3n^2+3n+1, -2n^3(n+1)] | Cor 2.2 | (0,0,1) |
| cat #4 | 24/(18G-11) = PCF[3n^2+3n+1, -2n^3(n+2)] | Cor 2.3 | (0,0,2) |
| cat #5 | 720/(450G-299) = PCF[3n^2+3n+1, -2n^3(n+3)] | Cor 2.4 | (0,0,3) |
| cat #6 | 2/(2G-1) = PCF[3n^2+7n+3, -2n^3(n+2)] | Cor 2.5 | (0,2,0) |
| cat #7 | 4/(2G+1) = PCF[3n^2+7n+3, -2n^2(n+1)(n+2)] | Cor 2.6 | (0,2,1) |
| cat #8 | 16/(6G-1) = PCF[3n^2+7n+3, -2n^2(n+2)^2] | Cor 2.7 | (0,2,2) |
| cat #9 | 288/(90G-31) = PCF[3n^2+7n+3, -2n^2(n+2)(n+3)] | Cor 2.8 | (0,2,3) |
| cat #10 | 1/(2-2G) = PCF[3n^2+9n+7, -2n(n+1)^3] | Cor 2.9 | Cohen family B (u,v,w)=(1,1,1); = tail of (0,0,γ→-1) |
| cat #11 | 24/(18G-11) = PCF[3n^2+11n+5, -2n^3(n+4)] | Cor 2.10 | (0,4,0) |
| cat #12 | 16/(6G-1) = PCF[3n^2+11n+5, -2n^2(n+1)(n+4)] | Cor 2.11 | (0,4,1) |
| cat #13 | 64/(18G+13) = PCF[3n^2+11n+5, -2n^2(n+2)(n+4)] | Cor 2.12 | (0,4,2) |
| cat #14 | 4/(6G-5) = PCF[3n^2+11n+9, -2n^2(n+2)^2] | Cor 2.13 | (2,2,0) |
| cat #15 | 8/(3-2G) = PCF[3n^2+11n+9, -2n(n+1)(n+2)^2] | Cor 2.14 | (2,2,1) |
| cat #16 | 32/(2G+5) = PCF[3n^2+11n+9, -2n(n+2)^3] | Cor 2.15 | (2,2,2) |
| cat #17 | 192/(18G+13) = PCF[3n^2+11n+9, -2n(n+2)^2(n+3)] | Cor 2.16 | (2,2,3) |
| cat #18 | 6/(17-18G) = PCF[3n^2+13n+13, -2n(n+1)^2(n+3)] | Cor 2.17 | Cohen family B (1,1,3) |
| cat #19 | 48/(90G-79) = PCF[3n^2+15n+15, -2n^2(n+2)(n+4)] | Cor 2.18 | (2,4,0) |
| cat #20 | 32/(19-18G) = PCF[3n^2+15n+15, -2n(n+1)(n+2)(n+4)] | Cor 2.19 | (2,4,1) |
| cat #21 | 128/(17-6G) = PCF[3n^2+15n+15, -2n(n+2)^2(n+4)] | Cor 2.20 | (2,4,2) |
| cat #22 | 8/(54G-49) = PCF[3n^2+15n+19, -2n(n+2)^3] | Cor 2.21 | Cohen family B (2,2,2) |
| cat #23 | 12/(83-90G) = PCF[3n^2+17n+23, -2n(n+1)(n+3)^2] | Cor 2.22 | Cohen family B (1,3,3) |

Status of all of #2–#23: **PROVED (Yamamoto 2024, arXiv preprint, unrefereed)**. Independent corroboration
for #2: Cohen DB 5.4.29 records it as obtained from/related to other labelled CFs by Apéry acceleration
("AP->5.4.31", "APD->5.3.1", "Dual to Pi") [F]. Other related work, not proofs: Naccache–Yifrach-Stav,
arXiv:2210.15669v5 — generalises the Catalan table to a two-parameter family Q_{c,κ} but states "we do not
claim them to be theorems" [F]; the three Catalan items in the Nature paper's Table 5 are cat #6, #14, #18
[F]. Contrary bookkeeping: site says "most of which are still unproven"; arXiv:2610.27799v1 Table 3 lists
"RM catalogue #1–#22 [88]" with the caption "Ramanujan Machine Catalan limits remain conjectural" [F].
CLAIMED on the forum: "yaron", May 2026, "1/2G" and "2 Proofs for Catalan's Constant G" — single page
images of longer documents (I read one: page 6 of 6, a 3F2 sketch for 2/(2G-1)) [F for that page];
unverifiable as posted.

### 1.3 pi^2 list — `.../2021/08/Zeta2_Results.pdf` [F]

12 entries: 1 "known", 2 "new and proven [1]", 9 "new and unproven". Site text: "most of which are
still unproven". All 9 are settled.

| id | value = PCF[a_n, b_n] (all num-checked) | site label | status today |
|---|---|---|---|
| pi2 #1 | 30/pi^2 = PCF[11n(n+1)+3, n^4] | known | Apéry |
| pi2 #2 | 8/pi^2 = PCF[3n(n+1)+1, -n^3(2n-1)] | new and proven [1] | PROVED: Kadyrov–Shynarbek–Orynbassar, arXiv:2103.03554 (v. of 2 Feb 2024), Thm 1.5 [F statement]; Yamamoto Cor 2.23 ((α,β,γ)=(0,0,-1/2)) [F]; C. Wang, arXiv:2602.03027v3 [A]; Cohen DB 5.4.6 [F] |
| pi2 #3 | 16/(4+pi^2) = PCF[3n(n+1)+1, -n^3(2n-3)] | new and unproven | **PROVED**: Yamamoto Cor 2.24, with a second Petkovšek-type proof at the end of §3 [F] ((0,0,-3/2)) |
| pi2 #4 | 24/pi^2 = PCF[7n(n+1)+2, 8n^4] | new and unproven | **PROVED (classical)**: this is the Franel-number Apéry limit. Franel numbers A(n)=sum C(n,k)^3 satisfy (n+1)^2 A(n+1) = (7n^2+7n+2)A(n) + 8n^2 A(n-1); with B(0)=0, B(1)=1, lim B/A = pi^2/24 was observed by Cusick (van der Poorten 1979) and **proved by Zagier 2009** ("Integral solutions of Apéry-like recurrence equations") — as stated in Chamberland–Straub, arXiv:2011.03400, §7 eq. (20) [F]. The bridge P_n=(n+1)!^2 A(n+1), Q_n=(n+1)!^2 B(n+1) is a one-line normalisation [I]. Cohen DB 5.4.10 [F]. Not in Yamamoto. Forum: "yaron" post "24/(pi^2)" (same vacuous template as his zeta(3) posts [I]). |
| pi2 #5 | 18/pi^2 = PCF[n(5n+6)+2, -2n^3(2n-1)] | new and proven [1] | PROVED: Kadyrov–Shynarbek–Orynbassar Thm 1.6 [F statement]; Cohen DB 5.4.7 [F] |
| pi2 #6 | 16/(pi^2-4) = PCF[n(3n+7)+3, -n^2(n+2)(2n-1)] (the site prints b_n with a shifted index, -(n+1)^2(n+3)(2n+1); the displayed fraction 3 - 3/(13 - 48/(29 - 225/...)) fixes the reading) | new and unproven | **PROVED**: Yamamoto Cor 2.25 ((0,2,-1/2)) [F]. CLAIMED also by Lezhe Gao, forum post 18 Apr 2026 (Bauer–Muir sketch, no details) [F] |
| pi2 #7 | 32/pi^2 = PCF[n(3n+7)+3, -n^2(n+2)(2n-3)] | new and unproven | **PROVED**: Yamamoto Cor 2.26 ((0,2,-3/2)) [F] |
| pi2 #8 | 16/(pi^2-8) = PCF[n(3n+11)+9, -n(n+2)^2(2n-1)] | new and unproven | **PROVED**: Yamamoto Cor 2.27 ((2,2,-1/2)) [F] |
| pi2 #9 | 16/(12-pi^2) = PCF[n(3n+11)+9, -n(n+2)^2(2n+1)] | new and unproven | **PROVED**: Yamamoto Cor 2.28 ((2,2,1/2)) [F] |
| pi2 #10 | 32/(32-3pi^2) = PCF[n(3n+15)+15, -n(n+2)(n+4)(2n+1)] | new and unproven | **PROVED**: Yamamoto Cor 2.29 ((2,4,1/2)) [F] |
| pi2 #11 | (16+3pi^2)/(16-pi^2) = PCF[n(3n+9)+7, -(n+1)^3(2n-3)] | new and unproven | **PROVED**: Yamamoto Cor 2.30 (tail of (0,0,-5/2): = 3/(PCF[3n^2+3n+1, -n^3(2n-5)] - 1)) [F] |
| pi2 #12 | 18/(pi^2-8) = PCF[n(5n+14)+10, -2n^3(2n+3)] | new and unproven | **PROVED**: Yamamoto Thm 3.1 (explicit solution y_n = 2^n n! (2n+3)!!) [F, proof read]; Cohen DB 5.4.7 sporadic entry (5n^2+4n+1, -2n^3(2n+3)) [F] |

Cohen DB parametric families containing these (stated experimentally) [F]: label 5.4.6
`a=3n^2+(8u+4v-3)n+(2u-1)(2u+2v-1), b=-n(n+2u)(n+2u+2v)(2n+2u-2k-1)`; label 5.4.7 (four families for
5n^2..., -2n(...)(2n±1)); label 5.4.10 (families for 7n^2..., 8n(...)).
Forum: "yaron" post "16/(4+pi^2)" (June 2026) — same template as his zeta(3) posts [I].

### 1.4 Mixed zeta orders — `.../2022/07/results_different_zeta_orders.pdf` [F]

5 entries, all "New and unproven"; site: "all of them still unproven". All 5 are **PROVED** — they are
Euler-type (telescoping) PCFs.

| id | value = PCF[a_n, b_n] (num-checked) | status today |
|---|---|---|
| mix #1 | -1/(zeta(4)+4 zeta(2)-8) = PCF[n^4+(n+1)^4+2(n^2+(n+1)^2), -n^8] | **PROVED**: Brier–Naccache–Yifrach-Stav, arXiv:2211.01058v2, §2 lemma (1/L = sum 1/(f(n+1)g(n)g(n+1))) + Example 2.1 [F]; Yamamoto Thm 3.3 (y_n=(2n+1) n!^4) [F]; Cohen DB 5.8.9 with series [F]; L. Gao forum post 31 May 2026 (complete elementary proof) [F] |
| mix #2 | 2/(2 zeta(5)+6 zeta(3)-9) = PCF[n^5+(n+1)^5+6(n^3+(n+1)^3), -n^10] | **PROVED**: Yamamoto Thm 3.4 (y_n=(3n^2+3n+1) n!^5) [F]; Cohen DB 5.8.10 [F]; Gao, Zenodo 10.5281/zenodo.20527727 [F abstract+Thm]; method of BNY [F] |
| mix #3 | **2/(2 zeta(5)-2 zeta(3)+1)** = PCF[n^5+(n+1)^5+6(n^3+(n+1)^3)-4(2n+1), -n^10] | **PROVED**: Yamamoto Thm 3.5 (y_n=(n^2+n+1) n!^5) [F]; Cohen DB 5.8.11 [F]; Gao Zenodo [F]. **Sign erratum**: the site PDF (and the statement of Yamamoto's Thm 3.5) print "-1" in the denominator; the num-check gives 2.98622... = 2/(2 zeta(5)-2 zeta(3)+1), which is also what Yamamoto's own displayed series sums to [I], what arXiv:2308.11829v2 Table 1 prints [F], and what Gao states [F] |
| mix #4 | 64/(64 zeta(5)+176 zeta(3)-273) = PCF[n^5+(n+1)^5+16(n^3+(n+1)^3)-4(2n+1), -n^10] | **PROVED**: Yamamoto Thm 3.6 [F]; Cohen DB 5.8.12 [F]; Gao forum comment 6 Jun 2026 [F] |
| mix #5 | 1/(zeta(7)-4 zeta(3)+4) = PCF[n^7+(n+1)^7+8(n^5+(n+1)^5)-8(n^3+(n+1)^3)+4(2n+1), -n^14] | **PROVED**: Yamamoto Thm 3.7 (y_n=(2n^2+2n+1) n!^7) [F]; Cohen DB 5.8.13 [F]; Gao forum comment [F] |

(The site PDF's displayed fractions also contain typos: it prints 1^8, 2^8, 3^8 for all rows and "275" for
a_1 = 75 in row 3 [F image]; the polynomial columns are right.) "yaron" posted six image "proofs" for
these in June 2026 — same vacuous template [I].

### 1.5 log 2 list — `.../2020/06/other.pdf` [F]

| id | value | status today |
|---|---|---|
| log2 #1 | 1/(1-log 2) = PCF[3n^2+7n+4, -2n^2(n+1)^2] (num-checked) | **PROVED**: Yamamoto Cor 2.31 ((α,β,γ)=(1,1,0)) [F]; L. Gao, forum post 28 Aug 2025 (complete elementary proof: P_n = 2^(n+1)(n+1)!(n+2)!) [F] and Zenodo 10.5281/zenodo.19696993 (16 May 2026) [A]; Cohen DB label 5.3.65, family `3n^2+(2(u+v+w)-3)n+..., -2n(n+u)(n+v)(n+w)` with (u,v,w)=(0,1,1) [F]. "yaron" also posted an image (May 2026) [A] |

Also in Yamamoto (not on a site list): Thm 3.2, 1/(2 log 2 - 1) = PCF[3n+3, -2n^2] [F].

### 1.6 Counts

| list | "unproven" on site | PROVED | CLAIMED only | OPEN |
|---|---|---|---|---|
| zeta(3) | 2 | 2 (Z1 refereed sources; Z2 preprint only) | 0 | 0 |
| Catalan | 22 | 22 (single preprint) | 0 | 0 |
| pi^2 | 9 (+2 already "proven") | 9 (8 preprint, 1 classical) | 0 | 0 |
| mixed zeta | 5 | 5 (several independent proofs) | 0 | 0 |
| log 2 | 1 | 1 | 0 | 0 |

### 1.7 Nature paper (arXiv:1907.00205v1–v4; Nature 590 (2021) 67–73)

v4 (30 Apr 2020) Appendix A, Tables 3–5 read [F]; v1–v3 grepped only [A]; Nature web page: abstract only [A].

* Table 5 (pi^2, zeta(3), Catalan): every row is on the site lists above → all settled.
* Table 4 (pi): no "unproven" rows in v4. v4 §5.1: "over the span of a few months, proofs for all the
  formulas in the original manuscript were presented"; App. F.2: Gauss CF proves "all the results we found
  so far for variations on pi" [F]. The general conjecture `1 + 1(2z-1)/(4 + 2(2z-3)/(7 + ...)) = 2^(2z+1)/(pi C(2z,z))`
  is proved in App. F.2.1 [F].
* Table 3 (e), three rows "new and unproven" (deg a = 2, deg b = 3; all num-checked):
  * e-i: 6e/(2e-3) = PCF[n(n+6)+7, -n^2(n+3)] — **no proof located**.
  * e-ii: 3/(4e-10) = PCF[(n+1)(n+3), +n^2(n+3)] (the table prints "-(n+3)n^2" but the displayed fraction and
    the num-check require "+") — in Cohen DB as label 5.10.5 with an explicit series for the reciprocal
    (sum 2(-1)^n/(n!(n^2+3n+1)(n^2+5n+5))) [F]; Cohen gives no proof, but a series form means it is
    Euler-type and routine [I].
  * e-iii: 4e/(2e-1) = PCF[n(n+3)+3, -n^2(n+2)] — **no proof located**.
  Searched: arXiv "Ramanujan Machine" listing, Lu 1907.05563, Kadyrov–Mashurov 1912.03214, Lynch 2012.12692,
  Dougherty-Bliss–Zeilberger 2004.00090 (their family is [an^2+bn+1 : -an^2-bn], deg 2/2), Yamamoto, Cohen DB
  grep for exp(1) with cubic b. Status: **OPEN as far as I could find** (low standing; after u_n = n! v_n
  the recurrence has degree-2 coefficients, i.e. confluent-hypergeometric territory [I]).
* Site PDFs "Previous results" (e: 13 + 12 conjectures; pi: 15 + 16 + 15), 2019 [F for two of the five,
  images; A for the rest]: computer-generated, redundant, all with deg a <= 1, deg b <= 2 (or period-2
  interlacings). Not tracked item by item. The classes are covered by Gauss/confluent CFs, Lu, Kadyrov–Mashurov,
  Dougherty-Bliss–Zeilberger (doubly infinite family [n+k : an]); the first item of `results_pi_0101`
  (-pi/4 = 1/(-1 + 1/(-4 + (-2)/(-7 + (-9)/(-10 + ...))))) is proved individually in C. Wang,
  arXiv:2601.08461v2 (and as a family in 2601.11892v3) [F abstracts].

### 1.8 Other sources of posted conjectures (not itemised exhaustively)

| source | what is posted | status |
|---|---|---|
| arXiv:2308.11829v2 (PNAS 121 (2024) e2321440121), §2.3 eqs (5)–(8), App. C Tables 2–4 | zeta(2) family n^2+(n+1)^2+α(α-1), -n^4; ZigZagZeta family eq. (8); polylog family; Table 4 sporadic rows | families are Euler-type and stated as provable/proved by the authors (App. C: "can be converted to infinite sums [31] using Euler's continued fraction formula"; polylog: "rigorously proven") [F]. Table 4 rows 1–4, 8 = mixed list (proved). **Table 4 rows 5–7 are unidentified → OPEN, see 2.1** |
| arXiv:2212.09470v2 (ICML 2023), App. A Tables 1–3 | a few dozen rows marked "new and unproven" (the word occurs 28 times in the text; rows not counted individually): regular-continued-fraction-pattern conjectures (signed/interlaced CFs with linear patterns) for e, e^2, sqrt(e), tan 1, Bessel ratios J_k(1)/J_l(1) | no proof literature located → OPEN as far as I could find [F tables skimmed; I for status]; Hurwitz-type patterns, low standing |
| arXiv:2412.16818v1 (NeurIPS 2024), §4.1 eqs (8)–(9), §4.3 | PCF[4, 4n^2+2n] = 2G_GA/(4G_GA-3) (first numerator 6), PCF[4, 4n^2+2n-2] = (4G_GA-1)/(3G_GA-2), a lemniscate-constant CF with b_n = 4n^2-2n; two infinite families (a_n = i+2j+1, b_n = n^2+(i+k)n → pi if k=1, log 2 if k=0; a_n=-2n+j-1, b_n=-n^2+jn+k → root of b_n) | "formal proofs for these formula hypotheses remain an open challenge" / "awaiting proof" [F]. No proof located. deg a <= 1, deg b = 2 → classical 2F1 territory [I]. Read the exact initial terms from the paper before using (my transcription of eqs (8)–(9) is from pdftotext and is not num-checked) |
| arXiv:2412.12361v2 (ICLR 2025), §5 and App. F | 75 "new" integer relations between C-transforms; Conjecture 2: C[n^2/(k^2(1-4n^2))] = (2/k)/(ln(k+1)-ln(k-1)) | the "New" zeta rows are of the mixed-zeta (Euler) type [I]; Conj. 2 is, up to an equivalence transformation, Gauss' CF for artanh(1/k) [I] — not an open problem in substance |
| arXiv:2502.17533v4 (NeurIPS 2025) | harvested pi formulas; Table 17 = clusters "not yet unified" with the pi-CMF | unification tasks, not truth conjectures [F] |
| site page "New Results by the Community" | image formulas (Dalbec–Gottlieb 2021 — later found in Lorentzen–Waadeland; Zhu He 2020; Michalski 2019; Berndt–Elimelech Hurwitz-zeta family 2022; Berndt–Gottlieb zeta(3) 2019; polylog family with linked proof PDF; Brier–Naccache–Yifrach-Stav family) | images not transcribed [A]; no status given on the page except the polylog proof |
| GitHub org RamanujanMachine (repos RamanujanMachine, ResearchTools, LIReC, euler2ai, Blind-Delta-Algorithm; shallow clones, 2026-11-24) | code, READMEs, `tests/conjectures_tests.py`, `results/massey.pdf` | **no file tracking conjecture/proof status** [F: grep over all md/txt/tex/json/csv]. The GitHub API was not reachable from this session; ASyMOB, ramanujantools etc. not inspected beyond README |
| LIReC / BOINC database ("thousands of new formulas", arXiv:2308.11829 §2.3) | not publicly itemised in a stable list I could reach | unknown; treat as unproved raw data [I] |

---

## 2. Ranked shortlist of what is open today

### 2.1 (highest standing, explicit) Three unidentified zeta(5)-level PCFs — arXiv:2308.11829v2 §3.6 Table 1 and App. C.3 Table 4 (rows 5–7); PNAS 2024

With sigma(n,d) := n^d + (n+1)^d and b_n = -64 n^10 in all three (value = a_0 - 64*1^10/(a_1 - 64*2^10/(a_2 - ...))):

* **T1**: a_n = 8 sigma(n,5) - 15 sigma(n,3) + 9(2n+1) → 1.20426... (a_0=2, a_1=156)
* **T2**: a_n = 8 sigma(n,5) - 12 sigma(n,3) + 7(2n+1) → 2.45174... (a_0=3, a_1=177)
* **T3**: a_n = 8 sigma(n,5) + 20 sigma(n,3) - 5(2n+1) → 22.8410... (a_0=23, a_1=429)

Paper, verbatim (§3.6): "Three of the discovered formulas in the table did not fit any combination of zetas
(in numerical tests using PSLQ) and are thus presented with their numerical values." And Table 4 caption:
"Lines 5 to 7 are only formula candidates as we did not find a closed-form expression for them, despite their
recursion formulas being strikingly similar to the others. We would expect them to yield a formula involving
constants from the same level of the hierarchy as zeta(5)." [F]
Task: identify the closed form and prove it. How checked: num-check reproduces the paper's digits
(T1 converges slowly: 1.20434 at depth 8000); not in Cohen DB v5 (grep for n^10 entries) [F]; not in the
Ramanujan Library App. F table [F]; no arXiv hit. Remark [I]: after dividing by 8 these are of the same
shape as mix #2–#4 with non-integer coefficients, so the Brier–Naccache–Yifrach-Stav / Petkovšek mechanism is
the natural first thing to test; if a polynomial g with a_n g(n) = 8[(n+1)^5 g(n+1) + n^5 g(n-1)] exists the
value is a telescoping sum, and half-integer roots would bring in log 2 / pi^k-type terms, which would explain
the failed zeta-only PSLQ. Not tested here.

### 2.2 General conjectures on conservative matrix fields (CMFs)

arXiv:2507.08138v2 (Weinbaum, Leibtag, Kalisch, Shalyt, Kaminer; v2 15 Jul 2026), end of the experimental
section (text lines following Fig. 5) [F]:

> "Conjecture 1. (Local constancy of the asymptotic flag). Given a CMF M and a direction v satisfying the
> hypotheses of Proposition 4.5, let V be an angular neighborhood of v on which the eigenvalues remain distinct
> in modulus. Then, for any u in V, B(x, u) = B(x, v)U, U is invertible upper-triangular."

> "Conjecture 2. (Continuous asymptotic height). Given a CMF M and a direction v satisfying the hypotheses of
> Theorem 5.1, the irrationality-measure sequence, in the sense of Definition 9, of the CMF ratio L^{p,q}_{x,v}(n)
> converges, and the resulting irrationality measure, denoted δ(v), is continuous as a function of v/|v|."

Related informal statements in arXiv:2308.11829v2 [F]: §1.3 "We put forward the conjecture that these matrix
fields have the potential to unify all polynomial continued fraction formulas associated with a specific
constant"; Fig. caption in §3: "We hypothesize that all continued fractions and infinite sums converging to a
certain constant can be derived from a single conservative matrix field of that constant"; §3.5–3.6: existence
of a CMF for zeta(5) is open ("we have not yet found such a matrix field"). For the zeta(3) CMF the
trajectory-independence of the limit is proved (David, arXiv:2303.09318v3, cited as [31]) [F for the PNAS
statement, A for David's proof]. Status: OPEN (no proof found; 2610.27799 uses CMFs but does not claim these).

### 2.3 Euler's constant — arXiv:2610.13919v1 (9 Oct 2026), §4

> "Conjecture 1 (Conjectured γ Apéry-like Approximations). For any r >= 3 and for any t in T_r (definition 20),
> theorem 3 and consequently theorem 1 hold."

Here T_r = {t(a,b) : k(b) >= 2}, t(a,b) = (a_1..a_{r-1}; b_1..b_r) with a in Z_{<=0}^{r-1}, b in N^r, k(b) the
multiplicity of min(b) (Def. 20); the authors prove only r = 3, t in T_3^* = {t(a,(1,1,1))} (Thm 1, Thm 3) [F].
Status: OPEN (one month old).

### 2.4 Factorial-reduction (FR) conjectures — arXiv:2111.04468v2 = Arnold Math. J. (2024), §2.3–§2.6 [F]

> "Conjecture 1.1. Given a polynomial b_n, there exists an a_n for which PCF[a_n, b_n] has FR if and only if
> b_n is decomposable to linear factors over Q."

> "Conjecture 1.2. For each b_n that is decomposable to linear factors over Q, there exists at least one
> infinite family of a_n s for which every PCF[a_n, b_n] has FR."

> "Conjecture 1.3. For any d in N, and B, x_1, ..., x_{2d} in Q, the following PCF[a_n, b_n] has FR:
> b_n = B * prod_{i=1}^{2d}(n - x_i), a_n = (B/m) * prod_{i=1}^{d}(n + 1 - x_i) - m * prod_{i=d+1}^{2d}(n - x_i)
> for any m in Q\{0, ±sqrt|B|}"

> "Conjecture 1.4. For every b_n of the form b_n = B(n - x_1)(n - x_2) with B, x_1, x_2 in Q, PCF[a_n, b_n] has
> FR if and only if a_n belongs to one (or more) of the families a_n^(k) = (B/m - m) n + k (Eq. 6),
> a_n^(k) = k(2n - x_1 - x_2 + 1) (Eq. 7), for any k in Q\{0}, and m in Q\{0, ±sqrt|B|}."

> "Conjecture 2 (The exact forms of GCD_n). We can represent every GCD_n as a multiple of two parts, a factorial
> and an exponential expression: ..." (long descriptive statement, §2.4–2.5)

> "Conjecture 3. For any PCF: GCD2_n ≐ GCD_n" (§2.6; ≐ means equal exponential order)

Also arXiv:2308.11829v2 §2: "We conjecture that the property of factorial reduction is a signature of
recurrence sequences that converge to mathematical constants." Status: OPEN, no proofs located. Remark [I]:
the family of Conj. 1.3 is exactly the Euler family a_n = h_1(n+1) + h_2(n), b_n = -h_1(n)h_2(n) (David,
arXiv:2308.02567), for which the convergents are explicit, so 1.3 looks accessible; the "only if" halves of
1.1 and 1.4 are the hard parts.

### 2.5 Convergence-rate conjecture — arXiv:2412.12361v2 §3.1

> "Conjecture 1. Let f_n be a real sequence. The formula C[f_n] converges in the following cases, at a rate
> provided by the error ε_n := |1 + f_1/(1 + f_2/(... + f_n))) - C[f_n]|. • If f_n = O(n^k) for some k < 0,
> then ε_n = O((n!)^k) ..." (further cases follow in the paper; quoted only as far as I read) [F partial]

The authors: "these discovered rules lack a formal proof for many of the special cases" (§6). OPEN.

### 2.6 Lower standing

* Nature Table 3: **e-i** 6e/(2e-3) = PCF[n(n+6)+7, -n^2(n+3)] and **e-iii** 4e/(2e-1) = PCF[n(n+3)+3, -n^2(n+2)]
  (no proof located; see 1.7).
* ESMA conjectures (arXiv:2212.09470 Tables 1–3), Gauss/lemniscate-constant PCFs and the two families of
  arXiv:2412.16818 §4.1/§4.3 (see 1.8).

### 2.7 Not open (do not spend effort)

All of zeta3 #4–#5, cat #2–#23, pi2 #3–#12, mix #1–#5, log2 #1; e-ii; Ramanujan Library Conjecture 2 [I].
If a refereed-quality proof matters to you, the items whose only proof is a preprint are: cat #2–#23,
pi2 #3, #6–#12, log2 #1 (Yamamoto — but mix-type and log2 also have independent elementary proofs) and
zeta3 #5 (Shvets).

---

## 3. Adjacent facts met during the survey

* arXiv:2610.27799v1 "New Directions in Irrationality" (RM group, 20 Oct 2026) claims mu(zeta(2)) < 4.62675,
  mu(zeta(3)) < 5.4197 and Q-linear independence of 1, zeta(3), zeta(4) [F abstract]; its introduction refers to
  a Zenodo preprint by A. Fauzan, "zeta(5) is irrational" (Sept 2026, doi 10.5281/zenodo.22826419, ref. [19]) and
  to an "OpenAI Math Release preprint, September 24, 2026" titled "Catalan's constant is irrational" (ref. [20])
  [F]. I did not read either preprint and take no position on them.
* `github.com/openai/math` (README, CONTENTS.md, history.md fetched 2026-11-24; "722 manuscripts organized
  into 374 families") [F README; CONTENTS.md grepped]: **no manuscript concerns Ramanujan Machine conjectures
  or polynomial continued fractions** (zero hits for "continued fraction"; the "Ramanujan" hits are
  Ramanujan conjecture / Ramanujan graphs). Adjacent families: 005 "Irrationality of Catalan's constant",
  017 "The irrationality exponent of pi is 2". The README says the collection is at "different stages of
  verification"; I did not read those manuscripts.
* Chamberland–Straub, arXiv:2011.03400: Conjecture 9 (Apéry limit zeta(2)/(d+1) for sums of d-th powers of
  binomials, d >= 3; proved only for d = 3 by Zagier) — not an RM conjecture; status after 2021 not checked.

---

## 4. Sources read

| source | date / version | label |
|---|---|---|
| ramanujanmachine.com: Results, zeta3.pdf, catalan.pdf, Zeta2_Results.pdf, results_different_zeta_orders.pdf, other.pdf (text + rendered images) | pages as of 2026-11-24; PDFs 2020-06, 2021-08, 2022-07 | F |
| ramanujanmachine.com: e/pi PDFs (5 files, 2019) | 2019 | F for 2 (images), A for 3 |
| ramanujanmachine.com: "Prove Our Conjectures" forum index + idea pages of Lezhe Gao (5) | posts Aug 2025 – Jun 2026 | F |
| forum images of "yaron" (Chun Wing Lai): 8/(7 zeta(3)), 12/(7 zeta(3)), Catalan p.6 | May–Jun 2026 | F for these 3; the other ~12 posts I (same titles/template) |
| ramanujanmachine.com: Publications, Community results, CMF page, Riddles, home | 2026-11-24 | F (text; images not transcribed) |
| Raayoni et al., arXiv:1907.00205 v4 (main text, App. A, F) | 30 Apr 2020 | F; v1–v3 A; Nature page A |
| Yamamoto, arXiv:2403.09729v1 | 13 Mar 2024 | F (all statements; proofs of Thm 2.1, 2.2, 3.1–3.7) |
| Kadyrov–Shynarbek–Orynbassar, arXiv:2103.03554 | version of 2 Feb 2024 | F (abstract, Thm 1.5/1.6 statements) |
| Shvets, arXiv:2604.06239v1 | 5 Apr 2026 | F (§1–2), A (§3–6) |
| Shvets, arXiv:2604.14219v2 | 24 Apr 2026 | F (abstract, intro) |
| Golyshev, arXiv:0908.1458 | 2009 | F (Thm 2.4, §3) |
| Golyshev–Kerr–Sasaki, arXiv:2009.14762; Batut–Olivier 1980; Almkvist–van Straten–Zudilin 2008; Zagier 2009 | — | A (cited through Shvets / Chamberland–Straub) |
| Chamberland–Straub, arXiv:2011.03400 | 2020 | F (§6–7) |
| Yang, arXiv:0709.1968 | 2007 | A (grep for Domb/level 6) |
| Cohen, arXiv:2304.11727 | 2023 | F (§1, §11–17) |
| Cohen, arXiv:2409.06086v5 + TeX source database (1925 CFs) | 16 Oct 2026 | F (intro, warnings; database grepped for the labels cited) |
| Cohen, arXiv:2401.17720 | 2024 | F (§5.3–5.5) |
| Brier–Naccache–Yifrach-Stav, arXiv:2211.01058v2 | 3 Nov 2022 | F (§1–2) |
| Naccache–Yifrach-Stav, arXiv:2210.15669v5 | 18 Nov 2022 | F (abstract, §2) |
| Dougherty-Bliss–Zeilberger, arXiv:2004.00090v3 | 26 May 2020 | F (§1–2) |
| Lu 1907.05563; Kadyrov–Mashurov 1912.03214v2; Lynch 2012.12692 | 2019–2020 | A |
| C. Wang, arXiv:2601.08461v2, 2601.11892v3, 2602.03027v3 | Mar–Apr 2026 | F (abstracts, theorem statements) |
| Gao, Zenodo 10.5281/zenodo.20527726 (v1 = .20527727) | 3 Jun 2026 | F (abstract, Thm 3.1 + proof, Example 1) |
| Gao, Zenodo 10.5281/zenodo.19696993 | 16 May 2026 | A (metadata) |
| Ben David et al., arXiv:2111.04468v2 | 2021 / AMJ 2024 | F (§2.3–2.6, §3) |
| Razon et al., arXiv:2212.09470v2 | 2022 | A/F (tables skimmed) |
| Elimelech et al., arXiv:2308.11829v2 | 2023 / PNAS 2024 | F (§1–4, App. C) |
| David, arXiv:2303.09318v3, 2308.02567v2 | 2023 | A (grep; Problem 17 read) |
| Shalyt et al., arXiv:2412.16818v1 | Dec 2024 | F (§4) |
| Beit-Halachmi–Kaminer, arXiv:2412.12361v2 | Dec 2024 | F (§3.1, §5–6, App. F table part) |
| Raz et al., arXiv:2502.17533v4 | 2025 | A/F (abstract, §5, App. J captions) |
| Weinbaum et al., arXiv:2507.08138v2 | 15 Jul 2026 | F (abstract, conjectures) |
| Kalisch et al., arXiv:2610.13919v1 | 9 Oct 2026 | F (§2, §4) |
| Leibtag et al., arXiv:2610.27799v1 | 20 Oct 2026 | F (abstract, Table 3 caption and RM rows, A.2.6, refs) |
| GitHub: RamanujanMachine/{RamanujanMachine, ResearchTools, LIReC, euler2ai, Blind-Delta-Algorithm} | shallow clones 2026-11-24 (last commits 2021-11, 2026-08, 2025-01, 2025-07, 2025-05) | F (READMEs, file lists, grep) |
| github.com/openai/math README, CONTENTS.md, history.md | 2026-11-24 | F / grep |
| arXiv API listings: "Ramanujan Machine", "polynomial continued fraction(s)", "conservative matrix field(s)", "Apery limit", authors Kaminer, Leibtag, Shvets, Wang, Yamamoto | 2026-11-24 | F |

Not done: arXiv:2607.09721 and 2610.40551 (Ramanujan Challenge) were ignored as instructed; MathSciNet/zbMATH
were not available, so "no journal reference found" for Yamamoto and Shvets rests on arXiv metadata and one
web search; the Nature Supplementary Information PDF was not fetched.
