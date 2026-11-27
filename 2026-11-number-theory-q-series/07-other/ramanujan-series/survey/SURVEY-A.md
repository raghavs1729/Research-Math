# SURVEY-A: hypergeometric Ramanujan-type series for 1/pi^2, 1/pi^3, 1/pi^4 and relatives

Status as of 2026-11-24. All numbered identities below were recomputed with `verify_A.py` (output in
`verify_A.out.txt`, same directory; 81 checks, minimum agreement 65 digits, run time about 3 minutes);
"digits" = digits of agreement reported there (mp.dps = 80). Statements that were recorded but not
recomputed are flagged as such where they occur.

Notation. (a)_n = a(a+1)...(a+n-1) (and Gamma(a+n)/Gamma(a) for non-integer n).
`R[s1,...,sd; z; P]` means  sum_{n>=0} z^n (s1)_n...(sd)_n / (1)_n^d * P(n).
zeta = Riemann zeta; L_D(s) = sum_n (D/n) n^-s (Kronecker symbol). phi5 = ((sqrt5-1)/2)^5 = (5 sqrt5 - 11)/2,
Phi = ((1+sqrt5)/2)^5 = 1/phi5.

Sources are cited by arXiv id where one exists. Short keys used throughout:

| key | reference |
|---|---|
| G02 | J. Guillera, Some binomial series obtained by the WZ-method, Adv. Appl. Math. 29 (2002); arXiv:math/0503345 |
| G03 | J. Guillera, About a new kind of Ramanujan-type series, Experiment. Math. 12 (2003) 507-510 (read from the EMIS mirror PDF) |
| G06 | J. Guillera, Generators of some Ramanujan formulas, Ramanujan J. 11 (2006); arXiv:1104.0392 |
| G08 | J. Guillera, Hypergeometric identities for 10 extended Ramanujan-type series, Ramanujan J. 15 (2008); arXiv:1104.0396 |
| G11 | J. Guillera, A new Ramanujan-like series for 1/pi^2, Ramanujan J. 26 (2011); arXiv:1003.1915 |
| G13a | J. Guillera, WZ-proofs of "divergent" Ramanujan-type series (2013); arXiv:1012.2681 |
| G13b | J. Guillera, More hypergeometric identities related to Ramanujan-type series, Ramanujan J. 32 (2013); arXiv:1104.1994 |
| G12 | J. Guillera, Kind of proofs of Ramanujan-like series; arXiv:1203.1255v3 |
| AG | G. Almkvist, J. Guillera, Ramanujan-like series for 1/pi^2 and string theory, Exp. Math. 21 (2012); arXiv:1009.5202 |
| GZ12 | J. Guillera, W. Zudilin, "Divergent" Ramanujan-type supercongruences, PAMS 140 (2012); arXiv:1004.4337 |
| G14s | J. Guillera, slides "Proofs of Ramanujan series by the WZ-method", Rutgers, 18 Sep 2014 (sites.math.rutgers.edu/~zeilberg/expmath/GuilleraSlides.pdf) |
| G15 | J. Guillera, A family of Ramanujan-Orr formulas for 1/pi; arXiv:1501.06413v4 |
| G16 | J. Guillera, Bilateral sums related to Ramanujan-like series; arXiv:1610.04839v3 |
| G18 | J. Guillera, Dougall's 5F4 sum and the WZ-algorithm, Ramanujan J. 46 (2018); arXiv:1611.04385v5 |
| G19 | J. Guillera, Bilateral Ramanujan-like series for 1/pi^k and their congruences, IJNT 16 (2020); arXiv:1908.05123 |
| G24 | J. Guillera, Heuristic derivation of Zudilin's supercongruences for rational Ramanujan series; arXiv:2312.16827v14 (May 2024) |
| G25 | J. Guillera, The WZ method and flawless WZ pairs; arXiv:2503.00570 (Mar 2025; cited by Au as Bull. Aust. Math. Soc.) |
| G25z | J. Guillera, On a formula for zeta(5) conjectured by Almkvist (PDF on anamat.unizar.es/jguillera/preprints, 2025) |
| CG | H. Cohen, J. Guillera, Rational hypergeometric Ramanujan identities for 1/pi^c: survey and generalizations; arXiv:2101.12592 |
| Z08 | W. Zudilin, Ramanujan-type formulae for 1/pi: a second wind?; arXiv:0712.1332 |
| DPVZ | Dembele, Panchishkin, Voight, Zudilin, Special hypergeometric motives and their L-functions: Asai recognition; arXiv:1906.07384v3 |
| Z26 | W. Zudilin, A supercongruence fantasy on Fibonacci and Lucas (after Guillera); arXiv:2602.06679 (Feb 2026) |
| Au23 | K. C. Au, Wilf-Zeilberger seeds and non-trivial hypergeometric identities; arXiv:2312.14051v4 (10 Oct 2024); J. Symbolic Comput. 130 (2025) per G25 |
| Au24q | K. C. Au, q-analogues of Wilf-Zeilberger seeds and Ramanujan 1/pi^k-formulas; arXiv:2403.04555v4 (29 Jul 2026) |
| Au22 | K. C. Au, Multiple zeta values, WZ-pairs and infinite sums computations; arXiv:2212.02986v4 |
| Au26 | K. C. Au, Discovering hypergeometric series with harmonic numbers via WZ seeds; arXiv:2602.08721v4 (9 Nov 2026) |
| Sh26a | N. Shalev, Borwein-Bradley-type identities for Ramanujan-type series for 1/pi^4; arXiv:2610.04770 (Oct 2026) |
| Sh26b | N. Shalev, A symmetric extension of Zhao's series for zeta(5); arXiv:2610.34482 (Oct 2026) |
| HS26 | Q.-H. Hou, Z.-W. Sun, Evaluations of some series via the WZ method; arXiv:2604.15172 (Apr 2026) |
| Sun26 | Z.-W. Sun, Various conjectural series identities; arXiv:2603.29973v3 (13 Apr 2026) |
| Qi26 | W.-W. Qi, Proof of two supercongruences of Guillera and Zudilin; arXiv:2604.04006 (Apr 2026) |
| FH25 | L.-Q. Feng, Q.-H. Hou, Finding congruences with the WZ method; arXiv:2506.19221 (Jun 2025) |
| MO | Y. Zhao, mathoverflow.net/questions/281009 (2017) |

Correction to the brief: arXiv:2209.01731 is an unrelated statistics paper. Au's WZ-seeds paper is
arXiv:2312.14051 (Au23). Also, CG does **not** mark a proved/unproved status per row; it only says (Sec. 5)
"some of the identities (and the three known identities for c >= 3) are still conjectural" (Jan 2021).
Per-identity status below was assembled from the other papers.

---------------------------------------------------------------------------------------------------

## 0. Executive summary

* Cohen-Guillera list 13 convergent rational series (10 for 1/pi^2, 1 for 1/pi^3, 2 for 1/pi^4) and
  8 "divergent" ones. As of today:
  * **proved**: CG #1,#2,#3,#4 (Guillera, WZ, 2002-2010), #6 (1930n^2+549n+45; Au23), #11 Gourevitch 1/pi^3 (Au23; simpler
    proof G25), #12 Zhao 1/pi^4 and #13 Cullen 1/pi^4 (Au23; G25; Sh26a), and the upside-down (dual) forms of all 8 divergent ones (#14-#21).
  * **still open as far as I could determine: exactly five rational series, all for 1/pi^2**:
    CG #5 (252n^2+63n+5, z=-1/48), #7 (5418n^2+693n+29, z=-1/80^3), #8 (532n^2+126n+9, z=(3/5)^6),
    #9 (1640n^2+278n+15, z=-1/2^10), #10 (1920n^2+304n+15, z=1/7^4);
    plus the Almkvist-Guillera **golden-ratio** series (irrational z) and its divergent/upside-down companions.
* I did **not** find 15 independent open "base" identities in this family: there are 6 base series
  (5 rational + 1 golden). The open list in Section 1 reaches 25+ entries only by including the conjectured
  companions of these 6 (shifted values at x = 1/2, Taylor/"extended" expansions, supercongruences, bilateral
  forms, divergent and upside-down golden companions). They are separate statements but all hang on the same 6 series.
* No proof, and no AI-assisted claim of a proof, of any of the 6 was found (searches of 2026-11-24, Section 5).
  The most recent explicit "still conjectural" statements I saw: Z26 (Feb 2026) for the golden series; Sun26
  (Apr 2026) for #7, #9, #10; Au24q v4 (Jul 2026) and Au23 for #10; Guillera's homepage for #7.
* WZ obstruction as stated by Au (Au23, Sec. 3.4): with his seeds, only Seed1, Seed4, Seed10 can produce a
  G(n,0) whose "exponent" z lies outside the group G = <n^n : n in Z> in Q^x. z = 1/48, 1/80^3, (3/5)^6, 1/7^4
  are outside G (my computation from his criterion: 3-adic valuation 1; 5-adic valuation 3; 5-adic valuation 6;
  7-adic valuation 4). z = -1/2^10 = -(4^-5) **is** in G, so #9 is not excluded by this criterion - yet no seed
  match has been reported for it.

---------------------------------------------------------------------------------------------------

## 1. STILL OPEN (as far as determined)

### 1.1 Rational convergent series for 1/pi^2 (the five unproved Guillera / Almkvist-Guillera series)

| id | statement | source | status / how checked | digits | WZ note |
|---|---|---|---|---|---|
| O1-G252 | R[1/2,1/3,2/3,1/4,3/4; -1/48; 252n^2+63n+5] = 48/pi^2 | G03 eq (2-3) (PSLQ, 2003); Z08 (5.5); CG Sec.5 #5; DPVZ Table 2.10 #9; G19 Table 3 #9 | OPEN. G12 Sec. 2.1 lists the only proved 1/pi^2 series as 120-, 820-, 20-, 74-series; G14s: "unproved"; G19: "most of them remain unproved"; G24 Ex. 7.2 (May 2024): "Using the PSLQ algorithm we conjecture"; absent from all examples of Au23 v4 and Au24q v4; arXiv/web searches found no proof. Caveat: Sun26 Remark 5.17 writes "= 48/pi^2 by [9]" with [9] = CG (a survey without proofs) - I read this as loose wording, not a proof. | 80 | z=-1/48 not in Au's group G. Same z as the proved 1/pi series R[1/2,1/4,3/4; -1/48; 28n+3] = 16/(sqrt3 pi) (G03 (2-6); G08 Identity 7 has a WZ pair for that one). |
| O2-G1640 | R[1/2,1/4,3/4,1/6,5/6; -1/2^10; 1640n^2+278n+15] = 256/(sqrt3 pi^2) | G03 (2-2); Z08 (5.4); CG #9; DPVZ #14; G19 Table 3 #14 | OPEN. G14s "remain unproved"; Sun26 Remark 5.21 (Apr 2026): "was conjectured to have the value ... (cf. [9])"; not in Au23/Au24q; no proof found. | 80 | z in G, so Au's exponent obstruction does not apply; no seed/shift reported. Au24q: "none of the 1/pi-formulas ... with s = 1/6 can currently be derived using WZ techniques" (said for 1/pi). |
| O3-G5418 | R[1/2,1/3,2/3,1/6,5/6; -1/80^3; 5418n^2+693n+29] = 128 sqrt5/pi^2 | G03 (2-4); Z08 (5.6); CG #7; DPVZ (2.12) | OPEN. DPVZ prints "?="; Guillera homepage (samples page, fetched 2026-11-24): "the fastest of all my formulas ... It has not been proved yet"; G24 Ex. 7.4 "we conjecture"; Sun26 Remark 5.20 "was conjectured to have the value"; no proof found. | 80 | z not in G (5^3). Same z=-1/80^3 and same 5418 as the Chudnovsky-type R[1/2,1/6,5/6; -1/80^3; 5418n+263] = 640 sqrt15/(3 pi), which Au23 names explicitly as not reachable by his seeds. |
| O4-G1920 | R[1/2,1/8,3/8,5/8,7/8; 1/7^4; 1920n^2+304n+15] = 56 sqrt7/pi^2 | G03 (2-5); Z08 (5.7); CG #10; AG intro; DPVZ (6.3) | OPEN. Au23 intro: "more developments are necessary if we want to use it to prove ..." this series, printed with "?="; Au23 Sec. 3.4: "do not succumb under this technique"; DPVZ "?="; Sun26 Remark 5.22 "was conjectured"; no proof found. | 80 | z not in G (7^4). Companion with the same z: P17 below (proved 2017 by a modular/Orr-type method, not WZ). |
| O5-AG532 | R[1/2,1/3,2/3,1/6,5/6; (3/5)^6; 532n^2+126n+9] = 375/(4 pi^2)  (equivalently sum (6n)!/n!^6 (532n^2+126n+9)/10^(6n) = 375/(4 pi^2)) | AG eq (1) and Table 3 (2010, via Calabi-Yau mirror-map search); CG #8; G15 (24) | OPEN. G14s "unproved"; G15 calls it "challenging"; not in Au23/Au24q; no proof found. Caveat as for O1: Sun26 Remark 5.18 writes "in view of [9]". | 80 | z not in G (5^6). |

### 1.2 Irrational (golden-ratio) series and its companions

| id | statement | source | status / how checked | digits | WZ note |
|---|---|---|---|---|---|
| O6-GOLD | R[1/2,1/2,1/2,1/3,2/3; (3 phi5)^3; (32-216 phi5)n^2+(18-162 phi5)n+(3-30 phi5)] = 3/pi^2. Equivalent (Z26 eq (8)): sum binom(2n,n)^4 binom(3n,n) [20(61-27 sqrt5)n^2+9(101-45 sqrt5)n+3(56-25 sqrt5)] ((sqrt5-1)/2)^(15n)/2^(6n) = 3/pi^2 | AG Table 3 (row "~5", k=8/3); G12 Sec. 3.1 ("unproved series"); DPVZ Table 2.10 #16 | OPEN. Z26 (Feb 2026) prints it with "?=" and builds new supercongruence conjectures on it. | 79 | z irrational, minimal polynomial z^2+36828z-729 (G12). Outside any rational seed/shift scheme as stated; G12: "the only known convergent series for 1/pi^2 with a value of z which is not rational". |
| O7-GOLDdiv | "divergent" conjugate: analytic continuation of R[1/2,1/2,1/2,1/3,2/3; z; (2408+216 phi5)n^2+(1800+162 phi5)n+(333+30 phi5)] at z = (-3/phi5)^3 = -27 Phi^3 "=" 36/pi^2 | G12 Sec. 3.1; G16 eq (14) / Appendix (38) | OPEN (conjectured by PSLQ). Verified here by the Barnes integral (1/2 pi i) int prod (s_i)_s/(1)_s^5 P(s) Gamma(1+s)Gamma(-s)(-z)^s ds. | 80 | - |
| O8-GOLDud | sum_{n>=1} (-phi5/3)^(3n) (1)_n^5/((1/2)_n^3 (1/3)_n (2/3)_n) [(2408+216 phi5)n^2-(1800+162 phi5)n+(333+30 phi5)]/n^5 = (1125/4) sqrt5 L_5(3) - 448 zeta(3) | G12 Sec. 3.1 ("we have guessed") | OPEN. (G12 prints the sum from n=0; it must start at n=1.) | 77 | upside-down (dual) of O7 |
| O9-GOLD1pi | sum_{n>=0} (-1/(3 Phi))^(3n) (1/2)_n^5/((1)_n^3 (1/6)_n (5/6)_n) [(2408+216 phi5)n^2+(608+54 phi5)n+(35+3 phi5)]/(6n+1) = 3 Phi^(3/2)/pi | G16 eq (18) (derived from the conjectured bilateral form of O7) | OPEN / status unclear (no proof found; its rational analogue P18 was proved in 2026). The right-hand side was garbled in the text extraction; the form above was fixed numerically. | 80 | - |

Related conjectures on these recorded but **not** re-verified here: Z26 eqs (9)-(10) (Dwork-type supercongruences
for the "rationalised" Fibonacci/Lucas sums S3, S4 built from O6 and its formal conjugate).

### 1.3 Conjectured companions of the open series O1-O5 ("extended"/shifted/bilateral/p-adic)

Shifted sum (CG Sec. 6): S(x) = sum_{n>=0} P(n+x) H_{n+x} / (sign(a)^n |a|^(n+x)), a = 1/z,
H_{n+x} = prod_s Gamma(s+n+x)/Gamma(s) / Gamma(1+n+x)^5. So S(0) is the series itself.

(a) Values at x = 1/2 (CG table in Sec. 6.3, column "pi^3 S(1/2)/(a0 D sqrt D)"; stated there without proof):

| id | statement | status | digits |
|---|---|---|---|
| O13-CG5h | for O1: S(1/2) = (48/pi^3) * 4 log(2^5/3^3) | OPEN (base series open) | 79 |
| O14-CG7h | for O3: S(1/2) = (128 sqrt5/pi^3) * log(2^74 5^5/3^54) | OPEN | 79 |
| O15-CG8h | for O5: S(1/2) = (375/(4 pi^3)) * 2 asin(164833/5^8) | OPEN | 80 |
| O16-CG9h | for O2: S(1/2) = (256/(sqrt3 pi^3)) * log(3^9/2^14) | OPEN | 80 |
| O17-CG10h | for O4: S(1/2) = (56 sqrt7/pi^3) * 2 asin(2241857/2^25) | OPEN | 80 |

(The superscripts 5^8, 2^25, 2^74, 3^54, 3^9, 2^14 were lost in the text extraction of CG and were re-established
numerically.)

(b) Taylor expansions at x = 0 (CG Sec. 6.3; equivalent, up to normalisation, to Sun26 Conjectures 5.17, 5.18,
5.20, 5.21, 5.22 of April 2026, which state sum_k g^(j)(k) for j <= 4 or 5; and to Guillera's 2024 homepage
conjecture that the x^(2m+1) coefficient is a rational multiple of L(chi, m+1)):

S(x) = a0 D [ sqrt(D)/pi^2 + 0 x - a2 D sqrt(D) x^2 + 0 x^3 + a4 D^3 L_D(2) x^4 - a5 a0 D^4 L_D(3) x^5 + O(x^6) ]

| id | series | D | a0 | a2 | a4 | a5 | status | digits (worst of x^0..x^5) |
|---|---|---|---|---|---|---|---|---|
| T1-G252 | O1 | 1 | 48 | 3/2 | 157/4 | 91/9 | OPEN | 85 |
| T2-G5418 | O3 | 5 | 128/5 | 3/2 | 887/32 | 21/8 | OPEN | 72 |
| T3-AG532 | O5 | 1 | 375/4 | 4/3 | 40 | 6944/1125 | OPEN | 74 |
| T4-G1640 | O2 | 12 | 32/9 | 7/24 | 757/1152 | 3/16 | OPEN | 73 |
| T5-G1920 | O4 | 28 | 1 | 1/7 | 31/336 | 1/21 | OPEN | 65 |

(The displayed normalisation of CG - where exactly the factors D sit - was ambiguous in the extracted text; the
form above is the one that matches numerically. Spot checks against Sun26: (5.34) sum g22'' = -448 sqrt7 and
sum g22'''' = 55552 sqrt7 pi^2 agree with T5; (5.23) sum g17'' = -192, sum g17''' = -480 pi i agree with T1 after
multiplying by Sun's factor e^(i pi x).)

(c) Supercongruences (p-adic companions). For p prime, p not dividing the denominators,
S_p = sum_{n=0}^{p-1} z^n (s1)_n...(s5)_n/(1)_n^5 P(n):

  Zudilin type:  S_p == P(0) (D/p) p^2 (mod p^5);
  Zhao type (CG Sec. 8): S_p == P(0) (D/p) p^2 + C5 L_D(4-p) p^5 (mod p^6).

| id | series | D | C5 | status | exact check (verify_A.py) |
|---|---|---|---|---|---|
| SC1-G252 | O1 | 1 | -91/9 | OPEN (conjecture) | both hold for all primes 7 <= p <= 43 |
| SC2-G1640 | O2 | 12 | -15/2 | OPEN | both hold, 7 <= p <= 43 |
| SC3-G5418 | O3 | 5 | -35/216 | OPEN | both hold, 7 <= p <= 43 |
| SC4-G1920 | O4 | 28 | 35/8 | OPEN | both hold, 11 <= p <= 43 |
| SC5-AG532 | O5 | 1 | 1953/50 | OPEN | both hold, 7 <= p <= 43 |

Sources: Zudilin-type pattern for 1/pi^2: Zudilin 2009 (JNT 129, cited in AG Sec. 4, G19, G24; not read directly);
AG Sec. 4 states SC5 mod p^5; G19 eq (13) and Tables 3-5 (general conjecture); CG Sec. 6.3/8 (C5 values);
G24 Sec. 5-6 (extension S(nu p) == S(nu)(chi/p)p^2 (mod p^5) for all nu >= 1, and the L_p-refinement).
Status: I found no proof of any Zudilin-type mod p^5 supercongruence for a *convergent* 1/pi^2 series. What I saw
proved: GZ12 main theorem (mod p^5 for the divergent 10n^2+6n+1 at z=-4), Qi26 (mod p^5 for the divergent
172n^2+75n+9 at z=-27/16), FH25 (20n^2+8n+1 only mod p^2, "partially confirm conjectures by Sun"). Au24q (Jul 2026)
calls the p-adic versions "most of which are widely open conjectures". This paragraph is a partial picture, not a
complete survey of the supercongruence literature.

(d) Bilateral forms (G19 Sec. 2-4; G16; G24 Sec. 2). Recorded, **not numerically verified here** (the left tail
needs analytic continuation of a series in 1/z):
  * B1 (for O4, G19 Example 4.2): with A(n,x) = prod (s_i)_{n+x}/(1)_{n+x}^5 P(n+x) z0^(n+x),
    e^{-i pi x} cos(pi x) prod_{s != 1/2} (cos pi x - cos pi s)/(1 - cos pi s) * sum_{n in Z} A(n,x)
    = (56 sqrt7/pi^2) [1 + sum_{k=1,2} alpha_k (cos 2 pi k x - 1) + beta_k sin 2 pi k x],
    alpha_1 = -79/2, alpha_2 = 23/2, beta_1 = (5/2) i, beta_2 = -(3/2) i.  "All the content of this section is conjectured."
  * G19 conjectures in general that all alpha_k, beta_k are rational (proved there only for rank 3, i.e. 1/pi),
    and conjectures the "B-supercongruences" (7) with the epsilon_0 values of its Tables 3-5
    (epsilon_0 = 1, -4, 1, 1, -7 for O1, O5, O3, O2, O4).

### 1.4 Status unclear (no proof found, but adjacent results make a proof plausible or possibly existing)

| id | statement | source | why unclear | digits |
|---|---|---|---|---|
| O11-RO10 | sum_{n>=0} 2^(-6n) (1/10)_n(3/10)_n(7/10)_n(9/10)_n/((1/2)_n (1)_n^3) (2100n^2+1160n+63)/(2n+1) = 200/pi | G15 last display (PSLQ; "there is no known Orr-type factorization from which we can derive it") | 1/pi relative with denominators 10; conjectural in arXiv:1501.06413v4 (Dec 2017). The journal successor (Guillera, More Ramanujan-Orr formulas for 1/pi, New Zealand J. Math. 47 (2017)) could not be accessed (HTTP 503), so I cannot say whether it was proved there. | 80 |
| O12-S35 | sum_{n>=0} (3/5)^(6n) (1/3)_n(2/3)_n(1/6)_n(5/6)_n/((1/2)_n (1)_n^3) (133n^2+79n+6)/(2n+1) = 625/(32 pi) | G14s; G15 eq (22) (PSLQ; "we believe that it should be possible to prove") | companion of O5 with the same z. Its sibling with z = 7^-4 (P17) was "recently established by Guillera" in the NZJM 2017 paper according to DPVZ Sec. 6; O12 is probably in the same paper but I could not read it. | 80 |
| O18-CG6h | for the 1930-series (P5): S(1/2) = (384/pi^3) * 2 log 2 | CG Sec. 6.3 #6 | base series proved by Au23 Ex. VII through an identity with a free parameter a; I did not check whether a = 1/2 there gives this value. | 80 |
| O19-CG12h | for Zhao's series (P8): S(1/2) = 9216 zeta(3)/pi^7 | CG Sec. 6.3 #12 | Au23 Ex. II and Sh26a give parametric extensions; not checked whether they cover the half-shift. | 80 |
| O20-CG13h | for Cullen's series (P7): S(1/2) = (2048/15)/pi^4 | CG Sec. 6.3 #13 | same remark (Au23 Ex. I, G25 Sec. 4.1, Sh26a). | 80 |
| O21-Gz5 | sum_{n>=1} 64^-n (1)_n^7/(1/2)_n^7 [ 4(84n^3-110n^2+48n-7)/n^8 + 56(7n^2-5n+1)(3n-1)(H_{2n}-H_n)/n^7 ] = 4 pi^4 log 2 + 186 zeta(5) | G25z (expansion conjectured by Almkvist 2009; Guillera "justifies" it heuristically via p-adic supercongruences) | no rigorous proof in G25z; Au26 treats harmonic extensions of this type (some proved, some marked "?="), not matched one-to-one here. | 80 |

Also in this class (recorded, not listed as numbered identities): the Taylor-expansion and bilateral statements
attached to the *proved* series (CG Sec. 6.3 rows #1-#4, #6, #11-#13; G24 Examples 2.1-2.4; G19 Examples 4.1, 4.3, 4.4;
Sun26 Conjectures 5.13-5.16, 5.19) - several are WZ-provable from the a-extended identities of G08/G11/Au23/Sh26a,
but I did not establish which ones have actually been proved.

---------------------------------------------------------------------------------------------------

## 2. PROVED (with reference) - do not re-prove

### 2.1 Convergent series for 1/pi^c

| id | statement | conjectured in | proved by | digits |
|---|---|---|---|---|
| P1-G20 | R[1/2 x5; -1/4; 20n^2+8n+1] = 8/pi^2 | - | Guillera, WZ (G03 says proved in G02/G06; G12: "in 2003"); WZ pairs in G08 Identity 8 | 80 |
| P2-G820 | R[1/2 x5; -1/2^10; 820n^2+180n+13] = 128/pi^2 | - | Guillera 2002 (G02; G06); pairs in G02 Sec. 3, G08 Identity 9 | 80 |
| P3-G120 | R[1/2,1/2,1/2,1/4,3/4; 1/16; 120n^2+34n+3] = 32/pi^2 | - | Guillera 2002 (G02; G06); pairs in G02 Sec. 3, G08 Identity 10 | 80 |
| P4-G74 | R[1/2,1/2,1/2,1/3,2/3; 27/64; 74n^2+27n+3] = 48/pi^2 | - | Guillera 2010 (G11 eq (17)); q-analogue is mock modular (Au24q Ex. XX) | 80 |
| P5-G1930 | R[1/2,1/3,2/3,1/6,5/6; -(3/4)^6; 1930n^2+549n+45] = 384/pi^2 | Guillera 2010 (AG Table 3; G14s, where it is misprinted "1936") | **Au23 Example VII** (Seed2); shorter flawless-pair proof G25 Sec. 4.2. Au24q: "a long-standing conjecture until it was proved recently". | 80 |
| P6-GOUR | R[1/2 x7; 1/64; 168n^3+76n^2+14n+1] = 32/pi^3 | B. Gourevitch 2002 (G03 eq (4-1)) | **Au23 Example XI** (Dougall7F6 + Dougall5F4); simpler pair from Seed9 in Au24q (Ex. XVII; "appendix" of an earlier version per Au23 Remark 4.3 and G25); G25 Sec. 4.5 gives the generalisation P19 | 80 |
| P7-CULLEN | R[1/2 x7,1/4,3/4; 1/2^12; 43680n^4+20632n^3+4340n^2+466n+21] = 2048/pi^4 | J. Cullen 2010 | **Au23 Example I** (Dougall7F6; reduces to P2); G25 Sec. 4.1; new proof Sh26a | 80 |
| P8-ZHAO | R[1/2 x5,1/3,2/3,1/4,3/4; -27/256; 4528n^4+3180n^3+972n^2+147n+9] = 768/pi^4 | Y. Zhao 2017 (MO) | **Au23 Example II** (Dougall7F6; reduces to P1); new proof Sh26a | 80 |
| P9-SUN12 | sum (4/27)^n (1/2)_n^7/((1/6)_n(5/6)_n(1)_n^5) (92n^3+54n^2+12n+1)/(6n+1) = 12/pi^2 | Z.-W. Sun | Au23 Example VIII (Seed7) | 80 |
| P10-SHALEV | sum_{n>=0} (-1)^n (2n)!^10 (5n)!/(n!^15 (10n)!) * 3(5532n^4+5464n^3+2173n^2+416n+32)/((10n+1)(10n+3)) = 240/pi^2 | new | Sh26b eq (1.3) (companion of Zhao's zeta(5) series; AI-assisted, see Sec. 5) | 80 |
| P17-RO7 | sum 7^(-4n) (1/8)_n(3/8)_n(5/8)_n(7/8)_n/((1/2)_n(1)_n^3) (1920n^2+1072n+55)/(2n+1) = 196 sqrt7/(3 pi) | G15 eq (21) (PSLQ) | Guillera, NZJM 47 (2017) eq (1.6), **as reported by DPVZ Sec. 6** ("recently established by Guillera"); I could not read the NZJM paper itself | 80 |
| P18-HS | sum (-1/27)^n (1/2)_n^5/((1)_n^3(1/6)_n(5/6)_n) (28n^2+10n+1)/(6n+1) = 3/pi  (= sum (28k^2+10k+1) binom(2k,k)^5/((6k+1)(-64)^k binom(3k,k) binom(6k,3k))) | G16 eq (17); Sun | HS26 (abstract: "we prove ... conjectured by the second author"), WZ method | 80 |
| P19-GGOUR | for all complex k: sum_{n>=0} (1/2+3k)_n^2 (1/2-k)_n^2 (1/2+k)_n^3/((1)_n^2 (1+2k)_n^2 (1+k)_n^3) * (1/2)_k (1/6)_k^2 (5/6)_k^2/(1)_k^5 * 64^-n (729/1024)^k * (6n+6k+1)(28n^2+8n+1+12k^2+56nk+8k) = 32/pi^3 | - | G25 Sec. 4.5 (WZ + Carlson). In the extracted text the exponents 2 look like subscripts 2n/2k; the reading above is the one that verifies (checked at k = 0, 0.3, 1). | 80 (k=0.3) |

### 2.2 "Divergent" series (analytic continuation; Barnes integral) and their upside-down (dual) series

CG #14-#21. D-ids are the divergent statements, U-ids the convergent duals
sum_{n>=1} z^-n (1)_n^{2c+1}/prod(s_i)_n * P(-n)/n^{2c+1} = A * L(c+1).

| CG # | D-id: analytic continuation "=" | U-id: dual identity | proved by | digits (D / U) |
|---|---|---|---|---|
| 14 | D14: R[1/2 x5; -4; 10n^2+6n+1] "=" 4/pi^2 | U14: sum_{n>=1} (-1/4)^n (1)_n^5/(1/2)_n^5 (10n^2-6n+1)/n^5 = -28 zeta(3) | D14: G13a eq (9) (WZ + Barnes). U14: equivalent to f(1/2) = 7 zeta(3) in G08 eq (13) | 87 / 80 |
| 15 | D15: R[1/2 x5; -2^10; 205n^2+160n+32] "=" 16/pi^2 | U15: sum_{n>=1} (-1/2^10)^n (1)_n^5/(1/2)_n^5 (205n^2-160n+32)/n^5 = -2 zeta(3) | U15: Amdeberhan-Zeilberger 1997 (cited in G03, Au24q, Sh26b); also G08 eq (16) at a=1/2. D15: WZ entry given in G18 Sec. 6 | 80 / 80 |
| 16 | D16: R[1/2,1/2,1/2,1/3,2/3; -27; 28n^2+18n+3] "=" 6/pi^2 | U16: sum_{n>=1} (-1/27)^n (1)_n^5/((1/2)_n^3(1/3)_n(2/3)_n) (28n^2-18n+3)/n^5 = -14 zeta(3) | U16: Chu, as cited in G18 Sec. 6 ("[2, p. 525, Example 14]": sum (-1/27)^n (1)_n^5/((3/2)_n^3(4/3)_n(5/3)_n) (28n^2+38n+13) = (21/2) zeta(3); the bibliographic entry is missing from the arXiv reference list) | 80 / 80 |
| 17 | D17: R[1/2,1/5,2/5,3/5,4/5; -5^5/2^8; 483n^2+245n+30] "=" 80/pi^2 | U17: sum_{n>=1} (-2^8/5^5)^n (1)_n^5/((1/2)_n(1/5)_n(2/5)_n(3/5)_n(4/5)_n) (483n^2-245n+30)/n^5 = -896 zeta(3) | U17: **Au23 Example IX** (Seed3). D17: G25 Sec. 5 | 88 / 80 |
| 18 | D18: R[1/2,1/3,2/3,1/4,3/4; -27/16; 172n^2+75n+9] "=" 48/pi^2 | U18: sum_{n>=1} (-16/27)^n (1)_n^5/((1/2)_n(1/3)_n(2/3)_n(1/4)_n(3/4)_n) (172n^2-75n+9)/n^5 = -1792 zeta(3) | U18: Chu-Zhang, Math. Comp. 83 (2014) Example 42, in the form sum (-16/27)^n (1)_n^5/((3/2)(4/3)(5/3)(5/4)(7/4))_n (172n^2+269n+106) = 63 zeta(3); WZ re-proof G18 Sec. 3. D18: WZ entry in G18 Sec. 6 | 80 / 80 |
| 19 | (z = 64, sign ambiguity; not verified) | U19: sum_{n>=1} 64^-n (1)_n^7/(1/2)_n^7 (21n^3-22n^2+8n-1)/n^7 = pi^4/8 | **Au23 Example XI** eq (4.2) (conjectured Guillera 2003, G03: H(1/2) = pi^4/2) | - / 80 |
| 20 | (z = 27/4; not verified) | U20: sum_{n>=1} (4/27)^n (1)_n^7/((1/2)_n^5(1/3)_n(2/3)_n) (92n^3-84n^2+27n-3)/n^7 = 8 pi^4 | **Au23 Example IV** (Seed7); G25 Sec. 4.3 | - / 80 |
| 21 | D21: R[1/2 x5,1/5,2/5,3/5,4/5; -5^5/2^10; 5532n^4+5600n^3+2275n^2+425n+30] "=" 1280/pi^4 | U21: sum_{n>=1} (-2^10/5^5)^n (1)_n^9/((1/2)_n^5(1/5)_n(2/5)_n(3/5)_n(4/5)_n) (5532n^4-5600n^3+2275n^2-425n+30)/n^9 = -380928 zeta(5) | U21 (Zhao, MO 2017): **Au23 Example X** (Seed7 + Example XV); four-parameter extension Sh26b. D21: G25 Sec. 5 | 87 / 80 |

Other proved relatives used as calibration:
* U22: sum_{n>=1} (16/27)^n (1)_n^3/((1/2)_n(1/3)_n(2/3)_n) (11n-3)/n^3 = 8 pi^2 (G13a eq (14); conj. Sun). 80 digits.
* P11-G18: sum_{n>=0} (27/64)^n (1)_n^3(1/6)_n(5/6)_n/(1/2)_n^5 (74n^2+101n+35)(6n+1)/(2n+1)^5 = 16 pi^2/3 (G11 eq (18)). 80 digits.
* P12-P15: half-shifts of P1-P4: S(1/2) = (8/pi^3)(14 zeta(3)/pi^2), (128/pi^3)(2 zeta(3)/pi^2), (48/pi^3)(2 pi/3),
  (32/pi^3)(pi/3) - equivalent to f(1/2) = 7 zeta(3), 256 zeta(3) (G08 eqs (13),(16)), G11 eq (18), and
  f(1/2) = 16 pi^2/3 (G08 eq (19)). 80 digits each.
* P16: half-shift of Gourevitch = 8/pi^3, equivalent to U19. 80 digits.

### 2.3 Exactly what Au proved and what he says is out of reach

Au23 (arXiv:2312.14051v4, 10 Oct 2024; journal: J. Symbolic Comput.), abstract: "we prove some hypergeometric
identities conjectures due to Z.W. Sun, J. Guillera and Y. Zhao etc., including two Ramanujan-1/pi^4, one 1/pi^3
formulas as well as a remarkable series for zeta(5)". By example:

| Au23 example | identity | seed stated in Au23 |
|---|---|---|
| I | Cullen 1/pi^4 (P7); reduces to P2 by the a^-5 Laurent coefficient | Dougall7F6, A=B=C=D=E=-1, K=2, a=...=e=1/2, k0=1 |
| II | Zhao 1/pi^4 (P8); reduces to P1 | Dougall7F6, A=1, B=C=D=E=0, K=1, a=...=e=1/2, k0=1 |
| IV | U20 (= 8 pi^4), conj. Guillera | Seed7 |
| VII | P5 (1930n^2+549n+45 -> 384/pi^2), conj. Guillera | Seed2 (g(n)=0, lim_n sum_k F = 0; uses sum (-1/8)^n (1/2)_n^3/(1)_n^3 (6n+1) = 2 sqrt2/pi) |
| VIII | P9 (Sun, 12/pi^2) | Seed7 |
| IX | U17 (= 896 zeta(3) up to sign), "conjecture appeared in [CG]" | Seed3 |
| X | U21 (Zhao's zeta(5)) | Seed7; needs Example XV (sum_k = 2(93 zeta(5) - 7 pi^2 zeta(3))) and U14 |
| XI | P6 (Gourevitch) and U19 (pi^4/8) | Dougall7F6 for F1; the g(n)-sum is G2(n,0) of a Dougall5F4 pair F2 |
| III, V, VI, XII-XVII | conjectures of Z.-W. Sun (harmonic-number twists, L-values) | Gauss2F1, Seed4, Dougall5F4, Seed1, Seed3, Seed10 |

What Au23 says remains open / out of reach:
* Intro: "this approach ... is still inchoate and more developments are necessary if we want to use it to prove more
  non-trivial identities like" the 7^-4 series O4 (printed with "?=").
* Sec. 3.4: "From the lists of seeds, only Seed1, Seed4, Seed10 can generate series whose exponent of G(n,0) is not
  in G [the subgroup of Q^x generated by n^n]. Therefore many (known or conjectural) series, for instances"
  R[1/2,1/6,5/6; -1/80^3; 5418n+263] = 640 sqrt15/(3 pi) and O4 "do not succumb under this technique."
* Footnote 1: P6, P7, P8 "are all the currently known 1/pi^3, 1/pi^4 formulas, according to [CG]".

Au24q (v4, 29 Jul 2026): q-analogues for P1-P4, P6, P7, P8 (21 examples). Says: "It does not yield WZ-style proofs,
and hence q-analogues, for all Ramanujan-type 1/pi-formulas. This includes the series (1.1) [26390n+1103], as well as"
R[1/2,1/4,3/4; -1/18^2; 260n+23] = 72/pi and R[1/2,1/4,3/4; 1/3^4; 10n+1] = 9/(2 sqrt2 pi); "none of the 1/pi-formulas
in (1.2) with s = 1/6 can currently be derived using WZ techniques. This limitation may indicate the existence of
further, as yet undiscovered, WZ-seeds." Mentions P5 only as proved; none of O1-O5 appears.
Au22 (Dec 2022): c = 0 identities (Sun's conjectures); says for c >= 1 the method "seems not yet mature enough".
Au26 (v4, Nov 2026): harmonic-number extensions; identities whose g(n) or limit term is not evaluated are marked "?=";
none of O1-O6 is addressed. arXiv:2610.17240 and arXiv:2509.19960 (Au, 2025-26) concern holonomic sequences from
meromorphic modular forms and L-values/elliptic-K moments; a grep for the O1-O5 polynomials found nothing relevant
(only abstract and keyword hits were read).

---------------------------------------------------------------------------------------------------

## 3. KNOWN WZ PAIRS / SEEDS in the literature for this family (calibration data)

Convention unless stated: F(n+1,k) - F(n,k) = G(n,k+1) - G(n,k). "validated" = checked numerically at several
integer points in `verify_A.py` Section W (ids W1-W12).

### 3.1 Guillera's explicit pairs (all validated)

W1 (G08 Identity 8, proves sum_n G(n+a,0) form of P1: f(a) = 8a sum (1/2)_n^4/(a+1)_n^4 (4n+2a+1)):
  B(n,k) = (-1)^n/(2^(12n) 2^(8k)) * (2k)!^4 (2n)!^5/((n+k)!^4 k!^4 n!^6)
  F = B * 8n(2n+4k+1),   G = B * (20n^2+8n+1+24kn+8k^2+4k)

W2 (G08 Identity 8, second pair):
  B(n,k) = (-1)^k (-1)^n/(2^(12n) 2^(6k)) * (2n)!^4 (2n+2k)! (2k)!^2 (n-k-1/2)!/(n!^7 (n+k)!^3 k!^2 (n-1/2)!)
  F = B * 32 n^3/(2n-2k-1),   G = B * (20n^2+12kn+8n+2k+1)

W3 (G08 Identity 9, for P2):
  B(n,k) = (-1)^n/(2^(20n) 2^(8k)) * (2n+2k)!^4 (2n)!^5/((2n+k)!^4 (n+k)!^4 n!^6)
  F = B * 128 n (6n+4k+1)
  G = B * [ (2n+2k+1)^4/(2n+k+1)^4 * (820n^2+180n+13+8k^2+20k+72nk)
            - (296nk^3+1056n^2k^2+1280n^3k+528n^4+800n^3+1344n^2k+608nk^2+28k^3+408n^2+384nk+40k^2+72n+16k+1) * 32nk/(2n+k+1)^4 ]

W4 (G08 Identity 9, second pair):
  B(n,k) = (-1)^k (-1)^n/(2^(20n) 2^(6k)) * (2n)!^4 (2n+2k)!^3 (n-k-1/2)!/(n!^7 (n+k)!^3 (2n+k)!^2 (n-1/2)!)
  F = B * 2048 n^3/(2n-2k-1)
  G = B * [ 820n^2+180n+13 + k/(2n+k+1)^2 * (1312n^3+1340n^2k+336nk^2+1456n^2+828nk+40k^2+472n+79k+36) ]

W5 (G08 Identity 10, for P3):
  B(n,k) = 1/(2^(16n) 2^(8k)) * (2k)!^2 (2n+2k)!^2 (4n)! (2n)!^2/((2n+k)!^2 (n+k)!^4 k!^2 n!^4)
  F = B * 32 n (4n+4k+1)
  G = B * [ 120n^2+34n+3+32k^2+128kn+16k + k (32n^3+8n^2k+16n+6kn+40n^2+k+2)/(2n+k+1)^2 ]
  WARNING: in the extracted text of the paper, G appears with extra denominators (4n+4k+1); that reading fails the
  WZ equation. The form above was reconstructed numerically (from F and G(n,0)) and satisfies it exactly.

W6 (G08 Identity 10, second pair):
  B(n,k) = (-1)^k/(2^(16n) 2^(6k)) * (2n)!^2 (4n)! (2k)!^3 (2n-k-1/2)!/(n!^6 (n+k)!^2 k!^3 (2n-1/2)!)
  F = B * 512 n^3/(4n-2k-1),   G = B * (120n^2+84kn+34n+10k+3)

W7 (G11 eqs (10)-(11)):
  U(n,k) = (1/2)_n^3 (1+k/2)_n (1/2+k/2)_n/((1)_n^3 (1+k)_n^2) * (1/2)_k^2/(1)_k^2
  F = U * 96 n^3/(2n+k),   G = U * 12k(8n^2+6kn+2n+k)/(2n+k);     sum_n G(n,k) = 48/pi^2 for k >= 1

W8 (G11 eqs (15)-(16), for P4; obtained from W7 by k -> k+n):
  U(n,k) = (27/64)^n (1/2)_n^3 (1/2+k)_n^2 (1+k/3)_n (1/3+k/3)_n (2/3+k/3)_n/((1)_n^3 (1+k)_n (1+k/2)_n^2 (1/2+k/2)_n^2) * (1/2)_k^2/(1)_k^2
  F = U * 96 n^3/(3n+k)
  G = U * [ n(2n+1)^2(74n^2+27n+3) + k P(n,k) ]/((n+k/3)(2n+k+1)^2),
  P(n,k) = (2n+1)(296n^3+164n^2+26n+1) + (480n^3+360n^2+78n+5)k + (176n^2+80n+8)k^2 + (24n+4)k^3;   sum_n G(n,k) = 48/pi^2

W9 (G02 Sec. 2; gives sum (-1)^n binom(4n,2n) binom(2n,n)^2 (20n+3)/2^(10n) = (1/2) sum binom(2n,n)^3 (42n+5)/2^(12n) = 8/pi):
  G = (-1)^n (-1)^k/(2^(10n) 2^(2k)) * (20n+2k+3) * binom(2k,k)^2 binom(2n,n)^2 binom(4n-2k,2n-k)/(binom(2n,k) binom(n+k,n))
  F = 64 (-1)^n (-1)^k/(2^(10n) 2^(2k)) * n^2/(4n-2k-1) * (same binomial quotient)

W10 (G02 Sec. 3; gives sum binom(4n,2n) binom(2n,n)^4 (120n^2+34n+3)/2^(16n) = (1/4) sum (-1)^n binom(2n,n)^5 (820n^2+180n+13)/2^(20n)):
  G = (-1)^k/(2^(16n) 2^(4k)) * (120n^2+84nk+34n+10k+3) * binom(2k,k)^3 binom(2n,n)^4 binom(4n-2k,2n-k)/(binom(2n,k) binom(n+k,n)^2)
  F = 512 (-1)^k/(2^(16n) 2^(4k)) * n^3/(4n-2k-1) * (same binomial quotient)
  Gamma-form extension to complex k (G02): G(n,k) = pi^-5 (120n^2+84nk+34n+10k+3) cos(pi k) Gamma(2n-k+1/2) Gamma(n+1/2)^6
  Gamma(k+1/2)^3/(Gamma(n+k+1)^2 Gamma(2n+1)^3)  [the tail of this display was cut in my extraction; taken as the analogue of the
  F display, which is complete: F = (512/pi^5) n^3/(4n-2k-1) * cos(pi k) Gamma(2n-k+1/2) Gamma(n+1/2)^6 Gamma(k+1/2)^3/(Gamma(n+k+1)^2 Gamma(2n+1)^3)].

W11 (G13a Example 1 eq (12); z = 4, dual gives U22):
  U(n,k) = (1/2)_n^2 (1+k/2)_n (1/2+k/2)_n 4^n/((1)_n^2 (1+k)_n^2) * (1/2)_k/(1)_k = (2n)!^2 (2n+k)! (2k)!/(n!^4 k! (n+k)!^2) /(16^n 4^k)
  F = U * 2n^2/(2n+k),   G = U * (6n^2+2n+k+4nk)/(2n+k)
  SIGN NOTE: as printed, these satisfy F(n+1,k)-F(n,k) = -(G(n,k+1)-G(n,k)).

W12 (G25 Sec. 4.4, "flawless" pair for sum_n G(n,k) = pi^2/8 for all k):
  U(n,k) = (1)_n (1+k)_n (1+2k)_n (1+k/2)_n (1/2+k/2)_n/((1/2)_n (1/2+k)_n (1+2k/3)_n (1/3+2k/3)_n (2/3+2k/3)_n) * (1)_k^2/(1/2)_k^2
  S(n,k) = n(2n+2k+3/2)/((1+2k)(2n+2k+1)(3n+2k+1)(3n+2k+2))
  R(n,k) = [22n^3+(34k+49)n^2+(12k^2+48k+35)n+8(k+1)^2]/(2(1+2n)(3+2k+3n)(2n+2k+1)(3n+2k+1)(3n+2k+2))
  F = U S z^n,  G = U R z^n.
  NOTES: the paper prints z = 4/27; the pair and the sum pi^2/8 only work with **z = 16/27** (checked: sum_n G(n,0) and
  sum_n G(n,3/2) equal pi^2/8 to 65 digits), and with the same opposite sign as W11: F(n+1,k)-F(n,k) = -(G(n,k+1)-G(n,k)).
  (Correspondingly Sec. 4.3 of G25 prints z = 16/27 for the pair whose k=0 series is U20 with z = 4/27; I did not test that pair
  because its R(n,k) is not printed.)

### 3.2 Guillera 2025 "flawless" data (F = U S z^n; R found by Zeilberger's Maple code; **not validated here**, R not printed)

* Cullen (G25 Sec. 4.1): U(n,k) = (1/2)_n^5 (1/2+k)_n^5 (1/2-k)_n (1+k/4)_n (1/4+k/4)_n (1/2+k/4)_n (3/4+k/4)_n/((1)_n^5 (1+k/2)_n^5 (1/2+k/2)_n^5) * (1/2)_k^4/(1)_k^4,
  z = 1/2^12, S(n,k) = 2^16 n^5 (4k+6n+1)(3+2k)(1+2k)/((4n+k)(4k^2+8k+3)(2k+1-2n)) (as printed). F(0,k) = F(+oo,k) = 0, so
  sum_n G(n,k) is 1-periodic; constant by Carlson; value 2048/pi^4 from k = 1/2 (only n = 0 survives). Replacing F(n,k) by F(n,k-n)
  "leads to" Zhao's series.
* 1930-series (G25 Sec. 4.2): U(n,k) = (1/2)_n^3 (1+k/3)_n (1/3+k/3)_n (2/3+k/3)_n (1/6+2k/3)_n (1/2+2k/3)_n (5/6+2k/3)_n/((1)_n^3 (1+k/2)_n^3 (1/2+k/2)_n^3) * (1/4)_k (3/4)_k/(1)_k^2,
  z = -3^6/4^6, S(n,k) = 768 n^3/(k+3n); constant sum 384/pi^2 obtained from k -> +oo.
* zeta(4) dual U20 (G25 Sec. 4.3): U(n,k) = (1)_n^3 (1+k)_n^3 (1+2k)_n (1+k/2)_n (1/2+k/2)_n/((1/2)_n^3 (1/2+k)_n^3 (1+2k/3)_n (1/3+2k/3)_n (2/3+2k/3)_n) * (1)_k^4/(1/2)_k^4,
  S(n,k) = n(2n+2k+3/2)/((1+2k)(2n+2k+1)^3(3n+2k+1)(3n+2k+2)); sum_n G(n,k) = pi^4/32.
* Gourevitch (G25 Sec. 4.5): pass from Au's pair to G~(n,k) = G(n-k,k); result is identity P19.
* Divergent D21 (G25 Sec. 5): F = U S z^n, z = -3125/1024,
  U(n,k) = (1/2)_n^5 (1/2+k)_n^4 (1/5+2k/5)_n (2/5+2k/5)_n (3/5+2k/5)_n (4/5+2k/5)_n (1+2k/5)_n/((1)_n^5 (1+2k)_n (1+k/2)_n^4 (1/2+k/2)_n^4) * (1/2)_k^4/(1)_k^4,
  S(n,k) = 640 n^5 (6n+4k+1)/((n+2k+1)(5n+2k))   (exponent placement read from a garbled layout; not validated).
* Guillera's Maple generators: anamat.unizar.es/jguillera/Maple.html (on-wz-pairs.txt, dougall-sum-wz.txt, WZ-RAMA-SEARCH-2024.mw).
  G18 gives Dougall-5F4 "entries" (a,b,c,d,e as linear functions of n) for: 205n^2+160n+32 at -2^10: `0+3n,0+2n,0+2n,0+2n,0+2n`;
  172n^2+75n+9 at -27/16: `1/2+3n,1/2+2n,1/2+n,1/2+n,0+3n`; Chu-Zhang Ex. 42: U(n,k) = A(3/2+3n, 1/2+0n, 1+n, 1+2n, 1+2n, k),
  F(n,k) = U(n,k)/(2n+1) * (1)_n^5/((3/2)_n(2/3)_n(4/3)_n(3/4)_n(5/4)_n) (-16/27)^n (as printed).

### 3.3 Au's seeds (Au23 Sec. 3.3; a WZ-seed f(a1,...,am,k) gives a WZ pair F(n,k) = f(a1+A1 n, ..., am+Am n, k+k0+K n))

Transcribed from the extracted text. "checked" = sum over k >= 0 equals 1 numerically (verify by hand; t12 in work/).

  Gauss2F1(a,b,c,k) = Gamma(c-a)Gamma(a+k)Gamma(c-b)Gamma(b+k)/(Gamma(a)Gamma(b)Gamma(k+1)Gamma(c+k)Gamma(c-a-b))   [also given as Mathematica code in Au23]
  Dougall5F4(a,b,c,d,k) = (a+2k)Gamma(a+k)Gamma(b+k)Gamma(c+k)Gamma(d+k)Gamma(a-b-c+1)Gamma(a-b-d+1)Gamma(a-c-d+1)
        / (Gamma(b)Gamma(c)Gamma(d)Gamma(k+1)Gamma(a-b+k+1)Gamma(a-c+k+1)Gamma(a-d+k+1)Gamma(a-b-c-d+1))        [checked]
  Dougall7F6(a,b,c,d,e,k) = (-1)^(a+e) (a+2k)Gamma(a+k)Gamma(b+k)Gamma(c+k)Gamma(d+k)Gamma(e+k)Gamma(a-b-c+1)Gamma(a-b-d+1)Gamma(a-c-d+1)
        Gamma(-a+b+c+e)Gamma(-a+b+d+e)Gamma(-a+c+d+e)Gamma(2a-b-c-d-e+k+1)
        / (Gamma(b)Gamma(c)Gamma(d)Gamma(e)Gamma(k+1)Gamma(-a+b+e)Gamma(a-b+k+1)Gamma(-a+c+e)Gamma(a-c+k+1)Gamma(-a+d+e)Gamma(a-d+k+1)
           Gamma(a-e+k+1)Gamma(a-b-c-d+1)Gamma(2a-b-c-d-e+1)Gamma(-a+b+c+d+e+k))
        [consistent, factor by factor, with Au's F for Examples I and II below]
  Seed2(a,c,d,k) = 2^(-2c-2k) (-1)^(a+c-d) Gamma(-a+2d-1)Gamma(a+2k)Gamma(d-c)Gamma(c+k)
        / (Gamma(a)Gamma(c)Gamma(k+1)Gamma(-a+d-1/2)Gamma(d+k)Gamma(-a-2c+2d-1)Gamma(a+c-d+k+3/2))      (a -> a/2, b -> a/2+1/2 in Balanced3F2)
  Seed3(a,b,d,k) = 2^(2d) (a+2k)Gamma(a-b+1/2)Gamma(a+k)Gamma(b+2k)Gamma(d+k)Gamma(2a-b-2d+1)
        / (Gamma(b)Gamma(d)Gamma(k+1)Gamma(a-b-d+1/2)Gamma(2a-b+2k+1)Gamma(a-d+k+1))                     [checked] (b -> b/2, c -> b/2+1/2 in Dougall5F4)
  Seed7(a,b,d,e,k) = (-1)^(a+e) (a+2k)Gamma(a-b+1/2)Gamma(a+k)Gamma(b+2k)Gamma(d+k)Gamma(e+k)Gamma(2a-b-2d+1)Gamma(-a+b+e+1/2)
        Gamma(-2a+b+2d+2e)Gamma(2a-b-d-e+k+1/2)
        / (Gamma(b)Gamma(d)Gamma(e)Gamma(k+1)Gamma(a-b-d+1/2)Gamma(-2a+b+2e)Gamma(2a-b+2k+1)Gamma(-a+d+e)Gamma(a-d+k+1)Gamma(a-e+k+1)
           Gamma(2a-b-d-e+1/2)Gamma(-a+b+d+e+k+1/2))                                                      (b -> b/2, c -> b/2+1/2 in Dougall7F6)
  Seed9(a,b,d,k) = (a+2k)(-1)^d 2^(4a-2b-2d) Gamma(a-b+1/2)Gamma(a+k)Gamma(b+2k)Gamma(d+2k)Gamma(-2a+2b+d+1)Gamma(-2a+b+2d+1)Gamma(2a-b-d+k)
        / (Gamma(b)Gamma(d)Gamma(k+1)Gamma(-a+d+1/2)Gamma(-2a+b+d)Gamma(2a-b+2k+1)Gamma(2a-d+2k+1)Gamma(-a+b+d+k+1))
        (b -> b/2, c -> b/2+1/2, d -> d/2, e -> d/2+1/2 in Dougall7F6)
  (Also in Au's list, not copied: Dixon3F2, Watson3F2, Balanced3F2, Seed1, Seed4, Seed10.)

Au's search recipe (Sec. 3.4): the "exponent" z and the Pochhammer part of G(n,0) equal those of F(n,0); scan A,...,K in [-N,N]
for the exponent, then (a,...,k0) in (1/l)Z/Z for the Pochhammer part; the polynomial R(n) "will always pop up automatically
when running Gosper's algorithm". Example given: Seed3(a+n, b-n, d+n, k+k0-n) has exponent 16/27. Central tool (Prop. 2.1):
sum_k F(0,k) + sum_n g(n) = sum_n G(n,0) + lim_n sum_k F(n,k), g(n) = lim_k G(n,k). Mathematica package:
researchgate.net/publication/376812577 (not downloaded).

### 3.4 Au's explicit F(n,k) (G(n,k) = F(n,k) * rational certificate found by Gosper; certificates mostly not printed)

  Ex. I (Cullen):  F = (2k+3n+1/2) Gamma(1-2n)^3 Gamma(n+1/2)^3 Gamma(k+n+1/2)^5 Gamma(k+4n)
                       / (Gamma(1/2-n)^7 Gamma(2n)^2 Gamma(k-n+3/2) Gamma(k+2n+1)^5)
                   [given verbatim as Mathematica code in Au23; equals Dougall7F6(1/2-n x5, k+2n) up to the sign factor]
                   G = F * (-2k+2n-1)(68728320 k^2 n^8 + 41932800 k n^9 + 11182080 n^10 + 59 more monomials)/(131072 n^5 (k+2n+1)^5 (4k+6n+1));
                   use n -> n+a and compare the a^-5 Laurent coefficients at a = 0.
  Ex. II (Zhao):   F = (-1)^n (1/2+2k+3n) Gamma(1-n)^3 Gamma(n+1/2)^3 Gamma(k+n+1/2)^4 Gamma(k+2n+1/2) Gamma(k+3n)
                       / (Gamma(k+3/2) Gamma(1/2-n)^3 Gamma(n) Gamma(2n) Gamma(k+n+1) Gamma(k+2n+1)^4)
                   [exponents fixed by matching Dougall7F6(1/2+n, 1/2, 1/2, 1/2, 1/2, k+n)]
  Ex. IV (U20):    F = n (2k+2n+3/2) Gamma(n+1)^3 Gamma(k+n+1)^3 Gamma(2k+n+1) Gamma(k+2n+1)
                       / (Gamma(k+3/2) Gamma(n+1/2)^3 Gamma(k+n+3/2)^3 Gamma(2k+3n+3))        (as printed; exponents read from layout)
  Ex. VII (P5):    F = 2^(-2k-6n) Gamma(-n-1) Gamma(n-1/2) Gamma(k+3n+2) Gamma(2k+3n+3/2)
                       / (Gamma(-n-5/2) Gamma(-n-1/2) Gamma(n) Gamma(n+1) Gamma(k+2n+1) Gamma(k+2n+2) Gamma(k+2n+3))   (as printed)
                   identity used: sum_k F(0,k+a) = sum_{n>=1} G(n-1,a), then the a^-1 coefficient; G(n,0) differs from the target by a telescoping term.
  Ex. VIII (P9):   F = (4k+4n+3) Gamma(n-1/2) Gamma(k+n+1/2) Gamma(k+n+3/2)^2 Gamma(2k+n+3/2) Gamma(k+2n+1)
                       / (Gamma(k+3/2) Gamma(-n-1/2)^2 Gamma(n)^2 Gamma(n+1) Gamma(k+n+1) Gamma(k+n+2)^2 Gamma(2k+3n+5/2))   (as printed; exponents read from layout)
  Ex. IX (U17):    F = (-4)^n (4k+6n+3) Gamma(n+1)^2 Gamma(2n+1) Gamma(2k+n+1) Gamma(k+2n+1)^2/(Gamma(n+1/2)^2 Gamma(k+n+3/2)^2 Gamma(2k+5n+3))
                   g(n) = 0, lim = 0; sum_k F(0,k+a) = sum_{n>=1} G(n-1,a); at a = 0 the k-sum is sum_{k>=0} (4k+3)/((2k+1)^3 (k+1)) (1)_k^2/(1/2)_k^2 = (7/2) zeta(3) (as printed).
  Ex. X (U21):     F = (-1)^n (4k+6n+3) Gamma(n+1)^4 Gamma(2n+1) Gamma(2k+n+1) Gamma(k+2n+1)^4/(Gamma(n+1/2)^6 Gamma(k+n+3/2)^4 Gamma(2k+5n+3))
                   g(n) != 0 (gives the 14 zeta(3) series U14).
  Ex. XI (P6):     F1 = (-1)^n (4k+6n+3) Gamma(n+1/2)^6 Gamma(n+1)^4 Gamma(k+n+1/2) Gamma(k+n+1)^4 Gamma(k+4n+1)
                        / (Gamma(2n) Gamma(2n+1/2)^4 Gamma(k-n+3/2) Gamma(k+2n+3/2)^4 Gamma(k+2n+2))             (as printed; exponents read from layout)
                   F2 = (3a+2k+3n+1) Gamma(a+n+1/2)^2 Gamma(a+n+1)^3 Gamma(a+k+n+1/2)^3 Gamma(a+k+n+1)
                        / (Gamma(-a-n+1/2) Gamma(2a+2n+1/2) Gamma(2a+k+2n+1) Gamma(2a+k+2n+3/2)^3)                (as printed; exponents read from layout)
                   The n-sum bisects into the series sum 2^(3-8n)(21n^3-22n^2+8n-1)(1)_n^6(1)_{2n}/(n^7 (1/2)_n^8); Laurent coefficient at a = 0 gives U19,
                   leading coefficient at a = -1/4 gives P6.
  Au24q Ex. XVII (Gourevitch, simpler): classical seed Seed9(1-2n, 1/2-n, 1/2-n, k+n) (q-version Seed9q2); set
                   F'(n,k) = F(n+k,k) + G(n+k,k+1), G'(n,k) = G(n+k,k); then F'(0,k) = 0 = lim_n F'(n,k), S(k) = sum_n G'(n,k) is constant,
                   evaluated at k = 1/2.
  Au24q Ex. XIII (Cullen): seed Dougall7F6(1/2-n x5, k+2n); F(0,k) = 0 = lim_n F(n,k); S(k) = sum_n G(n,k) constant; S(0) = G(0,1/2).
  Au24q closing remark: a q-analogue of P5 "could in principle be derived by considering Seed2q2(1/2-n, n, 1, 2n), but analytical
  subtleties seem even more pronounced".

Lines marked "as printed; exponents read from layout" come from a two-dimensional PDF layout in which exponent positions are
ambiguous; re-derive them from the seed + shift data before use (as was done here for Ex. I and II).

---------------------------------------------------------------------------------------------------

## 4. Discrepancies / typos found while verifying (all resolved numerically)

1. Brief: arXiv:2209.01731 is not Au's paper (it is arXiv:2312.14051).
2. Au23 Example IV: the displayed conjecture has numerator "92n^3-84n^2+21n-3"; the derivation in the same example and CG
   have 27n; numerics confirm **27n** (U20).
3. G18: Example 42 gives 63 zeta(3), Sec. 6 of the same paper quotes "69 zeta(3)"; numerics give **63**.
4. G25 Sec. 4.4 prints z = 4/27; the identity requires z = 16/27 (W12). Sign convention of G13a (12) and G25 Sec. 4.4 is
   opposite to the stated WZ equation.
5. G14s prints the 1930-series as "1936n^2+549n+45".
6. G12 Sec. 3.1 prints the upside-down golden series O8 as a sum from n = 0 (it starts at n = 1).
7. CG Sec. 6.3: the sum equals a0 D sqrt(D)/pi^2 (the sentence before the display says "a0 sqrt(D)/pi^2").
8. Sun26 (5.23) appears in my extraction as "sum g17'(k) = 0"; the CG expansion T1 times e^(i pi x) gives 48 i/pi for
   that first derivative (the second and third derivatives agree). I could not tell whether this is an extraction
   artefact or a misprint.

---------------------------------------------------------------------------------------------------

## 5. How the status was checked (all on 2026-11-24)

Papers downloaded from arxiv.org/pdf and read via pdftotext (in `dl/`): math/0503345, 0712.1332, 0904.0406, 1003.1915,
1004.4337, 1007.2290, 1009.5202, 1012.2681, 1104.0392, 1104.0396, 1104.1994, 1201.5233, 1203.1255, 1206.3981, 1501.06413,
1610.04839, 1611.04385, 1802.02023, 1803.08477, 1906.07384, 1908.05123, 1910.01961, 2001.08104, 2101.12592, 2212.02986,
2312.14051, 2312.16827, 2403.04555, 2404.16199, 2407.00621, 2502.15249, 2503.00570, 2504.01975, 2506.19221, 2509.05897,
2509.19960, 2602.04128, 2602.06679, 2602.08721, 2603.29973, 2604.03327, 2604.04006, 2604.15172, 2610.04770, 2610.17240,
2610.34482. Not every paper was read in full: the ones cited for a specific statement were read at that statement; the
others were grepped for the polynomials 252/5418/532/1640/1920/1930 and for "conjectur|open|unproved|proved".
Other sources read: G03 (EMIS mirror), G14s, Guillera's homepage pages publications/preprints/samples/Maple and G25z,
MathOverflow 281009.

arXiv API queries (export.arxiv.org/api/query, sorted by date, results inspected back to 2021-2022):
  au:"Kam Cheong Au"; au:"Jesus Guillera"; au:Zudilin; au:"Noam Shalev"; au:"John M. Campbell" AND (abs:Guillera OR abs:Ramanujan);
  au:"Chuanan Wei"; au:"Zhi-Wei Sun" AND (ti:series OR ti:pi); au:Zeilberger AND (abs:pi OR abs:Ramanujan OR abs:WZ);
  abs:Guillera; abs:Guillera AND (abs:conjecture OR ...); abs:supercongruence AND abs:Guillera; abs:"Ramanujan-like series";
  abs:"Ramanujan-type series"; abs:"Wilf-Zeilberger" AND (abs:pi OR abs:Ramanujan OR abs:zeta); abs:Gourevich OR abs:Gourevitch;
  abs:"hypergeometric motives" AND (...); abs:Asai AND abs:hypergeometric; abs:"Calabi-Yau" AND abs:Ramanujan AND abs:pi;
  ti:"1/\pi^2", abs:"1/\pi^3", abs:"1/\pi^4" (these TeX-string queries returned nothing - the API does not index them usefully).
Web searches (5): the polynomials "1920n^2+304n+15", "5418n^2+693n+29", "532n^2+126n+9", "1640n^2+278n+15", "252n^2+63n+5" with
"proof/proved/Asai/Hilbert modular/2025/2026"; "AI-assisted proof ... Guillera conjecture ... (GPT OR Claude OR AlphaProof OR
Ramanujan Machine OR Lean)". No proof of O1-O6 surfaced.

Limits of this check: (i) arXiv abstract search is keyword-based - a proof whose abstract mentions none of
Guillera/Ramanujan-type/WZ would be missed; (ii) no access to MathSciNet/zbMATH, journals behind paywalls, github.com, or
the NZJM 2017 paper; (iii) Zudilin's "Arithmetic hypergeometric series" (Russian Math. Surveys 2011), "More Ramanujan-type
formulae for 1/pi^2" (2007), Zudilin's 2009 supercongruence paper, Chan-Wan-Zudilin, Chu-Zhang 2014 and Chu's papers were
not read directly (only as cited by the papers above).

AI-assisted work (2024-2026) seen:
* Sh26a and Sh26b (N. Shalev, Oct 2026) state that "Large language models (Claude Opus 5.5 by Anthropic and GPT-6 by OpenAI)
  were used, under the author's direction" to explore Au's WZ seeds and work out proofs; certificates verified in Mathematica.
  They give parametric extensions and new proofs of the Cullen, Zhao and Zhao-zeta(5) series, prove ten harmonic-number
  conjectures of Sun for 1/pi^4 and pi^4 and five of six conjectures from Au26 Example X, and find P10. Sh26a "Further
  identities": "we have found numerically Borwein-Bradley-type identities for further Ramanujan-type series ... several open
  conjectures of Sun and others" - no claim about O1-O6.
* A web search result (not read) mentions a 2026 ICML paper by Shalyt et al. on a human-AI workflow unifying Ramanujan's 17
  series for 1/pi via conservative matrix fields; this concerns 1/pi, not 1/pi^2.
* No AI-assisted claim of a proof of any of O1-O6 was found.

---------------------------------------------------------------------------------------------------

## 6. Not covered

* Non-hypergeometric (Calabi-Yau / Ramanujan-Sato type) series for 1/pi^2 of AG Sec. 3.2 (Hadamard products A*alpha, B*epsilon,
  ..., and Almkvist's "transformations"), and Zudilin's quadratic-transformation series: not single-sum hypergeometric, not
  verified here, status not determined.
* Z.-W. Sun's large families of conjectures (Sun26 has over 150; series with harmonic numbers, trinomial coefficients, etc.):
  only the items that restate O1-O5 (Sun26 Conj. 5.17-5.22) or were proved by Au/Hou-Sun/Shalev are recorded.
* Quadratic-irrational 1/pi series (Aldawoud's thesis, cited in CG; Campbell-Cooper-Ye arXiv:2602.09352) - not 1/pi^c with c >= 2.
* Zudilin Z26 Fibonacci/Lucas supercongruences (5)-(11) and G19 B-supercongruences: recorded as conjectural, not re-verified.
* G03 lists other parameter families searched without success (e.g. (1/2,1/12,5/12,7/12,11/12), (1/2,1/10,3/10,7/10,9/10));
  no identities for them are known to me.
