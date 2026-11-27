# Zudilin's relations (6) for the 13 hypergeometric cases other than (1/2,1/2): what is proved (literature survey, 2026-11-26)

Scope. W. Zudilin, *A hypergeometric version of the modularity of rigid Calabi–Yau manifolds*, SIGMA 14 (2018) 086,
arXiv:1805.00544, eq. (6): `L(f,k)/F_k(1) ∈ Q` for k = 1, 2, 3, where `F_j(1)` are the ε-coefficients of the
Γ-normalised Frobenius series (his eq. (5)) at z = 1. The case (1/2,1/2) is settled in `THEOREM.md` and is not
treated again here. This file records, for the other 13 cases, (a) what is known about the newform f, and (b) whether
any of the three relations, or an equivalent statement about the conifold transition matrix, is **proved**.

Notation. `L_k := L(f,k)`. `ε_f` is the root number (sign of the functional equation). "BKSZ" = Bönisch, Klemm,
Scheidegger, Zagier, arXiv:2203.09426 (Comm. Math. Phys. 405 (2024) 134). "BGK" = Bönisch, Golyshev, Klemm, *Fibering
out Calabi–Yau motives*, arXiv:2510.03939. "LTYZ" = Long, Tu, Yui, Zudilin, arXiv:1705.01663 (Adv. Math. 393 (2021)
108058). "EHMM III" = Grove, Long, Rosen, Tu, arXiv:2607.25173 (July 2026). "RWZ" = Rogers, Wan, Zucker, arXiv:1303.2259
(Ramanujan J. 37 (2015)).

## 0. Summary table

The column "numerical constant" gives `L(f,k)/F_k(1)`. I recomputed these for this file: PARI/GP 2.15.4 `lfunmf` for
the L-values, with `lfuncheckfeq` ≤ −195, and my own evaluation of `F_j(1)` from Zudilin's eq. (5). The rational
identifications hold to 10^−28 or better. The constants are new in the sense that Zudilin and BKSZ do not print them.
BKSZ give the equivalent data `w^±/ω_f^±` in their Table 3.

| (r,t) | f: level, LMFDB, sign ε_f | k = 1 | k = 2 | k = 3 |
|---|---|---|---|---|
| (1/2,1/3) | 36.4.a.a, ε=+1 | numerical only (Zudilin; BKSZ Table 3); const −1/3 | numerical only; const 1/18 | numerical only; const −1/90. k=3 ⇔ k=1 by proved results (§1.3) |
| (1/2,1/4) | 16.4.a.a = η(4τ)^16/(η(2τ)^4η(8τ)^4), ε=+1 | **L-value proved** in closed form (RWZ (30)): `L_1 = (π/4)·4F3(½,½,½,½;1,1,1;1)`. Relation to F_1(1) **numerical only**; const −1/8 ⇔ `F_1^{(½,¼)}(1) = −2π·F_0^{(½,½)}(1)` | **L-value proved** (RWZ (31) + our THEOREM.md): `L_2 = −(π/2)S_1 = −(π/16)F_1^{(½,½)}(1)`. Relation numerical only; const 1/32 ⇔ `F_2^{(½,¼)}(1) = −2π·F_1^{(½,½)}(1)` | **L-value proved** (RWZ (30)): `L_3 = (π^3/32)·4F3(½^4;1^3;1)`. Relation numerical only; const −3/448 |
| (1/2,1/6) | 72.4.a.b, **ε=−1** | numerical only; const 1/3 | **trivially true**: `L_2 = 0` because ε_f = −1, and F_2(1) ≈ 60.54 ≠ 0. The meaningful statement is `F_2(1) ∈ Q·π√7·L(f⊗χ_{−7},1)` (numerical, this file) | numerical only; const −1/468 |
| (1/3,1/3) | 27.4.a.a, **ε=−1** | numerical only; const 1/9 | trivially true (L_2 = 0). Meaningful: `F_2(1) = −(π/2)·L(f⊗χ_{−4},1)` (numerical) | numerical only; const −1/243 |
| (1/3,1/4) | 9.4.a.a = η(3τ)^8 (CM), ε=+1 | **L-value proved** (Damerell-type; Li–Long–Tu SIGMA 14 (2018) 090 Cor. 5.4; EHMM III Thm A.2.1): `L_1 = Γ(1/3)^9/(32√3π^5)`. Relation numerical only; const −1/18 ⇔ `F_1(1) = −(9/16)Γ(1/3)^9/(√3π^5)` | **L-value proved**: `L_2 = Γ(1/3)^9/(96π^4)`. Relation numerical only; const 1/54 ⇔ `F_2(1) = (9/16)Γ(1/3)^9/π^4` | **L-value proved**: `L_3 = Γ(1/3)^9/(144√3π^3)`. Relation numerical only; const −1/216 ⇔ `F_3(1) = −(3/2)Γ(1/3)^9/(√3π^3)` |
| (1/3,1/6) | 108.4.a.a, **ε=−1** | numerical only; const 1/2 | trivially true (L_2 = 0). Meaningful: `F_2(1) ∈ Q·π√7 L(f⊗χ_{−7},1)` (numerical) | numerical only; const −1/504 |
| (1/4,1/4) | 32.4.a.a = η(4τ)^10/η(8τ)^2 − 8η(8τ)^10/η(4τ)^2, **ε=−1** | numerical only; const 1/8 | trivially true (L_2 = 0). Meaningful: `F_2(1) = −(2/3)·π√3·L(f⊗χ_{−3},1)` (numerical) | numerical only; const −3/1280 |
| (1/4,1/6) | 144.4.a.f, ε=+1 | numerical only; const −5/6 | numerical only; const 1/36 | numerical only; const −5/2304 |
| (1/6,1/6) | 216.4.a.c, **ε=−1** | numerical only; const 5/6 | trivially true (L_2 = 0); `F_2(1) ∈ Q·π√7 L(f⊗χ_{−7},1)` (numerical) | numerical only; const −5/4752 |
| (1/5,2/5) | 25.4.a.b, ε=+1 | **closest to proved**: BGK Thm 1 (2025) proves that the conifold entries w^± are integral-cocycle periods of `f_50 = 5f(τ) − 20f(2τ)`. With Manin's theorem this gives `w^+ ∈ Q^×·π^2 L_1`, but the rational number is not written down. A dictionary `F_1(1) ↔ w^+` (§1.4) is still needed. Numerical const −1/5 | same status, via w^−; numerical const 4/125 | same status as k=1; numerical const −3/625 |
| (1/8,3/8) | 128.4.a.b, **ε=−1** | numerical only; const 1/2 | trivially true (L_2 = 0); `F_2(1) = −(1/3)π√3 L(f⊗χ_{−3},1)` (numerical) | numerical only; const −3/2816 |
| (1/10,3/10) | 200.4.a.f, **ε=−1** | numerical only; const 7/10 | trivially true (L_2 = 0); `F_2(1) = −(10/33)π√3 L(f⊗χ_{−3},1)` (numerical) | numerical only; const −21/34000 |
| (1/12,5/12) | 864.4.a.a, **ε=−1** | numerical only; const −13/6 | trivially true (L_2 = 0); `F_2(1) = (6/673)π√7 L(f⊗χ_{−7},1)` (numerical) | numerical only; const −13/39744 |

For comparison, the settled case (1/2,1/2) (8.4.a.a) has constants −1/16, 1/32, −3/256, which agree with THEOREM.md.
That agreement checks the pipeline.

**Modularity (a).** It is proved for all 13 cases, with explicit forms:
* LTYZ Theorem 2 proves `H_p(α,β;1) = a_p(f_α) + χ_α(p)p` for p > 5, with f_α as in their Table 1 (labels above). The
  proof uses the modularity of rigid Calabi–Yau threefolds over Q (Gouvêa–Yui; Dieulefait–Manoharmayum), Faltings–Serre,
  and Beukers–Cohen–Mellit.
* LTYZ Theorem 1 proves the mod-p^3 supercongruences of Zudilin's Observation 1 for all 14 cases (Remark 4 in Zudilin).
* EHMM III Theorem 1.2 gives an independent explicit construction of f for the pairs `{1/d,1/e}`, d ∈ {2,3,4},
  e ∈ {2,3,4,6}, e | M_d.

**Relations (b).** For none of the 13 cases is any of the three relations (6) proved as a relation with `F_k(1)`.
There are three caveats.
1. **k = 2 is trivially true in 8 cases.** In (1/2,1/6), (1/3,1/3), (1/3,1/6), (1/4,1/4), (1/6,1/6), (1/8,3/8),
   (1/10,3/10) and (1/12,5/12), the root number is −1. So `L(f,2) = 0`, and Zudilin's `L(f,2)/F_2(1) ∈ Q` holds
   with ratio 0, provided `F_2(1) ≠ 0` (a certified numerical fact; all values are ≥ 31). This is a proof of the
   letter of (6), but the content is degenerate. The real period `F_2(1) ∝ w^− ∝ ω_f^−` is then a twisted value
   `π√|d|·L(f⊗χ_d,1)`, observed numerically above. Zudilin's (6) is literally correct but uninformative there.
   Neither Zudilin nor BKSZ remark on it: BKSZ only say "the periods of f can always be expressed by values at s = 1, 2
   of the associated L-function or its twists".
2. **k = 1 and k = 3 are equivalent in every case**, by proved results plus a finite computation:
   * Zudilin's Remark 5 states as rigorous that `F_3(1)/(π^2 F_1(1)) ∈ Q`. Numerically the constant is
     `csc^2(πr) + csc^2(πt) − 2/3`, to 10^−27 in all 14 cases; Zudilin does not state it.
   * Manin's periods theorem gives `L(f,3)/(π^2 L(f,1)) ∈ Q`, computable by modular symbols.
   * So the k = 1 and k = 3 relations stand or fall together.
3. For (1/2,1/4) and (1/3,1/4) the **L-values** themselves are proved in closed form (RWZ; CM/Chowla–Selberg). What is
   missing is the hypergeometric side, `F_k(1)`. The relations become explicit, purely hypergeometric identities
   (table). I found none of these identities in the literature.

## 1. Cross-cutting facts

### 1.1 Zudilin's own remarks (read in `scratchpad/p6/z.txt` and the rendered p. 4)
* Eq. (6) is introduced as "numerical calculations suggest conjectural inclusions".
* Remark 2: "None of the relations in (6) seem to be proved." Also: "when r = t = 1/2, we have an extra rational
  relation F_0(1) = F_2(1)/(2π^2) … It is the only case when F_0(1) is linearly dependent over Q with F_j(1)/π^j for
  j = 1, 2, 3."
* Remark 3: for (1/3,1/4), `L(η(3τ)^8,2) = Γ(1/3)^9/(96π^4)` and `L(η(3τ)^8,3) = Γ(1/3)^9/(144√3π^3)`, citing
  Damerell's strategy. This concerns the L-values only, not F_k.
* Remark 4: all 14 supercongruences are proved in LTYZ.
* Remark 5, the only rigorous relation among the F_j(1): "if z = 1 … both terms e^{−πimε} and e^{πimε} in the Fourier
  expansion of (9) do not show up. This allows one to rigorously establish that F_1(1) and F_3(1)/π^2 are rationally
  proportional – something that could follow from (6) complemented with the Manin–Shimura relation". No constant is
  given. My numerics give `F_3(1)/F_1(1) = π^2(csc^2πr + csc^2πt − 2/3)`. Its value is 4/3 for (1/2,1/2), 5/3 for
  (1/2,1/3), 7/3 for (1/2,1/4), 13/3 for (1/2,1/6), 2 for (1/3,1/3), 8/3 for (1/3,1/4), 14/3 for (1/3,1/6), 10/3 for
  (1/4,1/4), 16/3 for (1/4,1/6), 22/3 for (1/6,1/6), 10/3 for (1/5,2/5), 22/3 for (1/8,3/8), and 34/3 for (1/10,3/10).
* Table 1 of Zudilin gives eta-product expressions only for levels 8, 36, 16, 27, 9, 32, 144 and 25. Section 5 has
  weight-6 and Ramanujan-type analogues. There are no further results on the 14 cases.

### 1.2 BKSZ (arXiv:2203.09426; read: Intro, §2.5, §3.1, Table 2, Table 3, §4, Appendix B, Table 8)
* §3.1: the conifold transition matrix `T_µ` (3.6) has entries w^±, a^±, e^±, b, c, d. "Numerically we find that these
  are up [to] multiplication by rational numbers equal to w^±", and similarly e^± ↔ quasi-periods η_F^±.
* Table 3 lists `w^+/ω_f^+, w^−/ω_f^−, η_F^+/e^+, η_F^−/e^−` for all 14 levels: 8: 32, −256, 32, −4; 9: 108, −108, 36,
  −36; 16: 64, −256, 64, −16; 25: 100, −250, 100, −40; 27: 108, −486, 108, −24; 32: 256, −512, 256, −128; 36: 72,
  −432, 72, −12; 72: 432, −864, 432, −216; 108: 864, −1296, 864, −576; 128: 1024, −1024, 1024, −1024; 144: 1728,
  −1728, 1728, −1728; 200: 8000, −4000, 8000, −16000; 216: 5184, −2592, 5184, −10368; 864: 20736, −10368, 20736,
  −41472. All of these are numerical, "to 1000 digits" (footnote 18).
* What is proved: §4 builds a correspondence with a Kuga–Sato threefold **only for N = 8** (four quadrics, (1/2,1/2)):
  "provides a proof for the identification of the Calabi-Yau periods in this fiber with the periods and quasiperiods
  associated with the unique newform in S_4(Γ_0(8))". The full matrix (3.11) is "numerically" found even for N = 8.
* For the other cases the intro says only that "an intermediate result … can be obtained by the technique of
  fibering out of motives [11]", with [11] = BGK, "work in progress 2022" at the time.
* §2.5: `b = X^0(µ) = (2πi)^3 ϖ_0(µ)` is "in general not well understood". This b is the analogue of Zudilin's
  `F_0(1)`, which matches Zudilin's remark that F_0(1) is independent except for (1/2,1/2).

### 1.3 Proved reduction k = 3 ⇔ k = 1 (all cases)
This follows from two proved results:
* Remark 5 of Zudilin: F_3(1) ∈ Q π^2 F_1(1).
* Manin's periods theorem for weight 4: L(f,1), L(f,3)/π^2 lie in Q·(one period ω^+); the ratio is explicit by modular
  symbols.

So for each case it suffices to prove the k = 1 relation; k = 3 then follows after one finite computation. This is my
assembly; it is not stated in the literature.

### 1.4 Dictionary F_j(1) ↔ T_µ (needed to transfer any BKSZ/BGK-type result)
The F_j are Golyshev's Γ-normalised Frobenius solutions. Their values at z = µ are rational/π-power combinations of
the first column (0, σw^+ + w^−, b, w^+) of T_µ. The numerics (Zudilin plus my Table) say:
* F_1(1), F_3(1)/π^2 ∝ w^+;
* F_2(1) ∝ w^−;
* F_0(1) involves b.

A proof of this dictionary is a finite monodromy/Γ-class computation that I did not find written anywhere. For
(1/2,1/2) it was bypassed by direct analysis.

### 1.5 BGK (arXiv:2510.03939, Oct 2025; read in full, Appendix A rendered)
* Main theorem, **quintic (1/5,2/5) only**. Theorem 1: "The numbers w^±, e^± and a^± are periods of the modular forms
  f_50, F_50 and g_50", with `f_50 = 5f(τ) − 20f(2τ)` and f = 25.4.a.b. Precisely, ∫_∞^{γ∞}(2πi)^3(1,5τ,25τ^2)(f_50,
  F_50, g_50)dτ ≡ (r^+,r^−)(γ)·[[ω^+,η^+,α^+],[ω^−,η^−,α^−]] + boundary terms, where r^± are cocycles Γ_0^*(50) → Z^3.
  The quantities ω^± coincide with w^±, η^± = −e^±/4 and α^± = −a^±/10 (their eq. after (3)).
* Theorem 2: b, c, d are integrals of f_50, g_50, F_50 between the CM points τ_± = ±2/5 + i/10.
* Consequence (my reading): w^+ and w^− are proved to be rational multiples of the Eichler–Shimura periods of f, i.e.
  `w^+ ∈ Q^× π^2 L(f,1)` and `w^− ∈ Q^× π L(f,2)`. The rational numbers are computable from r^±(γ) for one γ but are
  not printed; BKSZ's 100, −250 are numerical.
* Appendix A lists for every family 3F2-to-4F3 integral identities (from Doran–Malmendier) and the level N of the
  cusp form obtained by pullback, "—" when the pullback is not a cusp form. It states "one proceeds as for our main
  example", and the abstract says the appendix lists "other families … for which the fibering out proves the
  modularity of the associated mixed period matrix". The levels are:
  * (1/2,1/2): 16, 16;
  * (1/3,1/4): 48, 48, 18;
  * (1/2,1/4): 64, 64, 48, 8;
  * (1/5,2/5): 75, 50;
  * (1/3,1/3): 27;
  * (1/4,1/4): 32, 32;
  * (1/2,1/3): —, 12;
  * (1/2,1/6): —, —, 72, —;
  * (1/3,1/6): 108, —, —;
  * (1/8,3/8): —, 128, 128, —;
  * (1/4,1/6): 72, —, —;
  * (1/10,3/10): 200, —;
  * (1/6,1/6): — (no cusp form);
  * (1/12,5/12): — (no cusp form).

  No details or theorems are given for these other families. So (1/6,1/6) and (1/12,5/12) are not even claimed to be
  reachable by this method.
* BGK cite K. Bönisch's Bonn PhD thesis (2023), "Modularity of special motives of rank four associated with Calabi–Yau
  threefolds", for the method. The Bonn repository blocked automated access, so I could not read it. The repository
  abstract (via search) says it treats four Calabi–Yau threefolds and contains "experimental identities". It should be
  checked by hand.

### 1.6 Explicit Hypergeometric-Modularity Method (Allen–Grove–Long–Tu I, II; Grove–Long–Rosen–Tu III)
* EHMM I (arXiv:2404.00711, eqs (1.2)–(1.3)) re-proves `(π^2/16)4F3(½^4;1^3;1) = L(f_8,…)`; this is level 8 only.
* EHMM II (arXiv:2411.15116) proves L-value formulas only for **weight-3** forms and 4F3(−1)/5F4(−1) values. It has no
  weight-4 rigid-CY L-values.
* EHMM III (arXiv:2607.25173, Jul 2026):
  * Thm 1.2 covers modularity only.
  * §6, "Future problems", Example 6.1 re-derives `(π/4)·4F3(½^4;1^3;1) = L(f16.4.a.a,1)` (= RWZ (30)) by a contour
    argument, and writes `F(κ(d,2),1)` as an Euler integral. Here `F(κ(4,2),1) = 4F3(¼,¾,½,½;1,1,1;1)` is
    Zudilin's F_0(1) for (1/2,1/4), i.e. a b-type value, not one of the F_k in (6).
  * Example 6.2 quotes a preprint in preparation (Hsu–Tu–Yang) expressing `L(f12.4.a.a,2)` by two 4F3(1)'s with
    lower parameters (1, 5/6, 5/6). Here f12.4.a.a is the χ_{−3}-twist partner of the (1/2,1/3) form 36.4.a.a; that
    L-value is not an L-value of 36.4.a.a.
  * They also write `L(η(3τ)^8,1)` as a combination of 4F3(1)'s with parameters (½,½,½,⅓;1,1,5/6) etc. These are
    non-Frobenius 4F3's, not F_k(1) for (1/3,1/4).
  * Appendix A.2, Thm A.2.1: `L(η(3τ)^8,3) = Ω/162 = (2π^2/9)L(η(3τ)^8,1) = (2√3π/9)L(η(3τ)^8,2)`, Ω = Γ(1/3)^6/Γ(2/3)^3,
    proved. This matches Zudilin's Remark 3.
  * Their open problems ask "How to use the partition (1.3) to obtain relations between periods of modular forms and
    hypergeometric values" and cite Zudilin's paper as an open direction.
* Barman–Maity (arXiv:2604.02723, 2026): explicit eta-quotient constructions, finite-field identities only, no
  L-values.

### 1.7 Other sources checked
* **Li–Long–Tu**, *A Whipple 7F6 formula revisited* (arXiv:2103.08858, Matematica 2022). Theorems 4–5 give modularity
  of the 4-dim pieces for pairs (c,f) including (1/2,1/3) [36.4.a.a] and (1/2,1/6) [72.4.a.b]. The archimedean
  Theorem 6 is **only for (1/2,1/2)** (6F5 and 7F6 as periods of f8.6.a.a, f8.4.a.a). Text: "We believe the appearance
  of g_{c,f} in the (1/2,1/3) case can be explained in a similar way. This is left to any interested reader." Whipple's
  formula relates 4F3(½,½,f,1−f;1,1,1;1) (F_0-type) to 7F6, so it does not touch F_k(1), k ≥ 1.
* **Li–Long–Tu**, SIGMA 14 (2018) 090 (arXiv:1803.06072). §4–5 (Cor. 5.4) prove `L(η(3τ)^8,3)` as a power of a beta
  value: the CM L-values for (1/3,1/4).
* **RWZ** (arXiv:1303.2259, §5, read). For `f_1 = η(4τ)^16/(η(2τ)^4η(8τ)^4)`, which is **the (1/2,1/4) form 16.4.a.a**
  (LTYZ Table 1, Zudilin Table 1):
  * (30): `L(f_1,3) = (π/2)L(f_2,2) = (π^2/8)L(f_1,1) = (1/8)∫_0^1 K(k)^2/√(1−k^2)dk = (π^3/32)4F3(½^4;1^3;1)`;
  * (31): `L(f_2,3) = (π/4)L(f_1,2) = (π^2/4)L(f_2,1) = … = (1/4)∫_0^1K^2 = (π^4/128)7F6(…)`.

  Our THEOREM.md gives `L(f_2,1) = −S_1/2`. Hence `L(f_1,1) = πS_0/4`, `L(f_1,2) = −πS_1/2` and `L(f_1,3) = π^3S_0/32`,
  where `S_0 = F_0^{(½,½)}(1)` and `S_1 = F_1^{(½,½)}(1)/8`. I checked all three against PARI to 40 digits. The
  reason: 16.4.a.a is the χ_{−4}-twist of 8.4.a.a (compare a_3: −4 vs +4).
* **D. Bao** (arXiv:2608.26514, §4) gives a period-polynomial explanation of the RWZ f_1/f_2 relations (Prop. 4.1:
  `F|σ = −4iG`, `G|σ = −(i/4)F`, σ = [[1,0],[4,1]]). There is nothing on the hypergeometric side. §5 concerns levels
  6/18, not our forms.
* **K. C. Au** (arXiv:2509.19960): dimension counts of L-values on Γ_1(4), Γ(4), Γ_1(8). Level-16 forms "are not
  twist-minimal, so can be ignored". Nothing on F_k(1).
* **Osburn–Straub** (arXiv:1806.05207): interpolated Apéry-like sequences ↔ L(f,2) for weight 3 and for level 8
  weight 4. Nothing for our 13 levels.
* **Golyshev–Zagier**, *Interpolated Apéry numbers, quasiperiods of modular forms, and motivic gamma functions*
  (PSPUM 103.2, 2021): level 8 only (the (1/2,1/2) case).
* **Dembélé–Panchishkin–Voight–Zudilin** (arXiv:1906.07384): 5F4 motives and Asai L-functions (Conjecture 5.1);
  not our 4F3 conifold fibres.
* **Candelas–de la Ossa–McGovern** (arXiv:2410.07107), CdlO–Elmi–van Straten (arXiv:1912.06146),
  Candelas–Kuusela–McGovern (arXiv:2104.02718): attractor points. 2410.07107 merely summarises BKSZ: conifold
  relations "observed (and even demonstrated rigorously in one case using a correspondence)". W. Yang
  (arXiv:2001.07211): rank-2 attractors of the CdlOEvS (AESZ 34) family; not conifold fibres of the 14 families.
* Rogers–Zudilin / Zudilin *Period(d)ness of L-values* (MPIM 12-42a) and Samart's Mahler-measure papers: weight 2
  (conductor 32 etc.) and weight-3/4 forms of levels 8/16 only, as far as the passages located show.
* Searches for papers citing Zudilin's SIGMA paper and for "conifold" + "L-value" + "proof" turned up nothing beyond
  the above.

## 2. Per-case notes

Throughout:
* "Relation" means `L(f,k)/F_k(1) ∈ Q`.
* "Numerical (Z)" means Zudilin 2018 asserts it numerically without constants. "Numerical (BKSZ)" means the
  equivalent statement `w^± ∈ Q ω_f^±` is in BKSZ Table 3 (numerical, 1000 digits).
* Constants are from this file's computation.
* The k = 3 relation is always equivalent to k = 1 (§1.3).

**(1/2,1/3)** — f = 36.4.a.a = η_6^14/(η_2^3η_18^3) − 3η_2^3η_6^2η_18^3 (Zudilin; LTYZ). Also f = χ_{−3}⊗f12.4.a.a
(EHMM III Thm 1.2, LLT 7F6 Thm 4). ε = +1. The L-values of 36.4.a.a are not known in hypergeometric closed form. The
Hsu–Tu–Yang result quoted in EHMM III concerns `L(f12.4.a.a,2)`, the untwisted partner, which is unpublished and not
the same L-function. BGK Appendix A: one cusp-form route (level 12), not carried out. Constants −1/3, 1/18, −1/90.
Status: k = 1, 2, 3 numerical only (Z; BKSZ).

**(1/2,1/4)** — f = 16.4.a.a = η_4^16/(η_2^4η_8^4) (Zudilin, LTYZ, RWZ's f_1). Also f = χ_{−4}⊗8.4.a.a and
f = η(2τ)η(4τ)^12/(η(2τ)^24 + 64η(4τ)^24)·(E_2(2τ) − 2E_2(4τ))^2 (EHMM III (1.4)). ε = +1.
* All three L-values are **proved** in terms of the (1/2,1/2) constants (RWZ (30), (31) + THEOREM.md):
  `L_1 = πS_0/4`, `L_2 = −πS_1/2`, `L_3 = π^3 S_0/32`.
* Zudilin's relations are therefore equivalent to the purely hypergeometric identities (numerically to 10^−27):

      F_1^{(½,¼)}(1) = −2π·F_0^{(½,½)}(1),     F_2^{(½,¼)}(1) = −2π·F_1^{(½,½)}(1),     F_3^{(½,¼)}(1) = −(14/3)π^3·F_0^{(½,½)}(1)

  The third identity follows from the first by Remark 5 once the constant 7/3 is pinned down.
* Constants −1/8, 1/32, −3/448.
* None of these transfer identities is in the literature. EHMM III §6 uses exactly the relevant contour/Euler-integral
  representation (`F(κ(4,2),1) = (2/π)∫_0^{1/2}(z/(1−z))^{1/2} 2F1(…;z)^2 dz/z`), but only for the ε^0 coefficient.
* BGK Appendix A lists four cusp-form routes (levels 64, 64, 48, 8), not carried out.

**(1/2,1/6)** — f = 72.4.a.b. LTYZ Table 1 gives a 4-term eta-quotient combination:
η_4^3η_6^6η_18^2/(η_12^2η_36) − 3η_2^2η_6^6η_36^3/(η_4η_12^2) + 8η_2^3η_12^6η_36^2/(η_6^2η_18) − 16η_12^12/η_6^4.
* ε = −1, so k = 2 is trivially true (ratio 0).
* F_2(1) = 60.5445…; numerically F_2(1)/(π√7 L(f⊗χ_{−7},1)) = 3/88 and F_2(1)/(π√11 L(f⊗χ_{−11},1)) = −3/224.
* Constants 1/3 (k=1), −1/468 (k=3).
* BGK route: one cusp form (level 72), not carried out. LLT 7F6: modularity only.

**(1/3,1/3)** — f = 27.4.a.a = η_1^3η_3^4η_9 − 27η_3η_9^4η_27^3 (Zudilin, LTYZ); also χ_{−3}⊗f27.4.a.b (EHMM III).
* ε = −1, so k = 2 is trivially true.
* Numerically `F_2(1) = −(π/2)L(f⊗χ_{−4},1)`; also = −(1/13)π√7 L(f⊗χ_{−7},1) = (1/30)π√8 L(f⊗χ_{−8},1).
* Constants 1/9, −1/243.
* BGK route: level 27, not carried out.

**(1/3,1/4)** — f = 9.4.a.a = η(3τ)^8, CM by Q(√−3). ε = +1.
* All three L-values are **proved**: Damerell/Chowla–Selberg; explicit in Zudilin Remark 3 (stated as known), Li–Long–Tu
  SIGMA 2018 Cor. 5.4, and EHMM III Thm A.2.1 (with proof).
* `L_1 = Γ(1/3)^9/(32√3π^5)`, `L_2 = Γ(1/3)^9/(96π^4)`, `L_3 = Γ(1/3)^9/(144√3π^3)`.
* Relations ⇔ closed Γ-evaluations of the ε-derivatives of 4F3(⅓,⅔,¼,¾;1,1,1;1) (numerically to 10^−27):

      F_1(1) = −(9/16)Γ(1/3)^9/(√3π^5),   F_2(1) = (9/16)Γ(1/3)^9/π^4,   F_3(1) = −(3/2)Γ(1/3)^9/(√3π^3)

  These are not in the literature.
* Constants −1/18, 1/54, −1/216.
* BGK routes: levels 48, 48, 18 (not carried out).

**(1/3,1/6)** — f = 108.4.a.a = η_6^10/η_18^2 − 27η_18^10/η_6^2 + 9η_6^7η_54^3/η_18^2 − 9η_2^3η_18^7/η_6^2 (LTYZ Table 1;
the 4th term as printed there, worth re-checking). Also χ_{−3}⊗f108.4.a.d (EHMM III).
* ε = −1, so k = 2 is trivially true.
* Numerically F_2(1)/(π√7 L(f⊗χ_{−7},1)) = −2/69.
* Constants 1/2, −1/504. BGK route: level 108.

**(1/4,1/4)** — f = 32.4.a.a = η_4^10/η_8^2 − 8η_8^10/η_4^2 (Zudilin, LTYZ). It is not CM (a_3 = −8, a_5 = −10).
Also χ_{−4}⊗f32.4.a.c (EHMM III).
* ε = −1, so k = 2 is trivially true.
* Numerically `F_2(1) = −(2/3)π√3 L(f⊗χ_{−3},1) = (1/20)π√7 L(f⊗χ_{−7},1)`.
* Constants 1/8, −3/1280. BGK routes: 32, 32.

**(1/4,1/6)** — f = 144.4.a.f = η_12^32/(η_6^12η_24^12) + 16η_6^4η_24^4 (Zudilin, LTYZ); also χ_{−4}⊗f72.4.a.d (EHMM III).
ε = +1. Constants −5/6, 1/36, −5/2304. All numerical only. BGK route: level 72.

**(1/6,1/6)** — f = 216.4.a.c; no eta-quotient formula found.
* ε = −1, so k = 2 is trivially true.
* Numerically F_2(1)/(π√7 L(f⊗χ_{−7},1)) = −6/241.
* Constants 5/6, −5/4752.
* BGK Appendix A: no cusp-form route ("—").

**(1/5,2/5)** (quintic) — f = 25.4.a.b = η_5^10/(η_1η_25) + 5η_1^2η_5^4η_25^2 (Zudilin, LTYZ, BKSZ (1.2)).
* ε = +1. Modularity of the conifold fibre: Schoen 1986.
* **BGK 2025 Thm 1 proves** that the conifold periods w^±, e^±, a^± are (cocycle-)periods of f_50 = 5f − 20f(2τ) and
  the associated meromorphic forms. The rational constants (BKSZ Table 3: w^+ = 100ω_f^+, w^− = −250ω_f^−) are not
  extracted in BGK, but are a finite computation from their monodromy data.
* To reach Zudilin's (6), one still needs (i) the dictionary §1.4 (F_1(1) ∝ w^+, F_2(1) ∝ w^−), and (ii) the Manin
  normalisation of ω_f^± against L(f,1), L(f,2).
* Constants −1/5, 4/125, −3/625.
* Status: "follows from proved results (BGK Thm 1 + Manin + Zudilin Remark 5) by a finite but unwritten argument".
  This is the only case besides (1/2,1/2) where the conifold period ↔ modular period identification is proved.

**(1/8,3/8)** — f = 128.4.a.b. LTYZ Table 1 gives a long 13-term eta-quotient combination.
* ε = −1, so k = 2 is trivially true.
* Numerically `F_2(1) = −(1/3)π√3 L(f⊗χ_{−3},1)`.
* Constants 1/2, −3/2816. BGK routes: 128, 128.

**(1/10,3/10)** — f = 200.4.a.f; no eta formula found.
* ε = −1, so k = 2 is trivially true.
* Numerically `F_2(1) = −(10/33)π√3 L(f⊗χ_{−3},1)`.
* Constants 7/10, −21/34000. BGK route: level 200.

**(1/12,5/12)** — f = 864.4.a.a; no eta formula found.
* PARI gives `L(f,2) = 0` to 10^−58, so ε = −1. This conflicts with the naive product of the Atkin–Lehner signs in
  BKSZ Table 8 (W_2 = W_3 = −1), which should be re-checked; the vanishing itself is solid (checkfeq −204).
* k = 2 is trivially true.
* Numerically F_2(1)/(π√7 L(f⊗χ_{−7},1)) = 6/673.
* Constants −13/6, −13/39744.
* BGK: no cusp-form route.

## 3. Judgement

* **Genuinely open (k = 1, equivalently k = 3), with no proved partial results beyond modularity:** (1/2,1/3),
  (1/2,1/6), (1/3,1/3), (1/3,1/6), (1/4,1/4), (1/4,1/6), (1/8,3/8), (1/10,3/10).
  * BGK claim their method applies to these, but nothing is written.
  * (1/6,1/6) and (1/12,5/12) are the hardest: there is no eta formula and no BGK cusp-form route.
* **k = 2 is (trivially) proved** for the 8 sign-(−1) cases. The meaningful replacement `F_2(1) ∈ Q π√|d| L(f⊗χ_d,1)`
  is numerical only (this file).
* **Closest to provable:**
  1. **(1/2,1/4)**. The L-values are proved (RWZ). The task reduces to two hypergeometric transfer identities:
     `F_1^{(½,¼)}(1) = −2πF_0^{(½,½)}(1)` and `F_2^{(½,¼)}(1) = −2πF_1^{(½,½)}(1)`, then k = 3 via the Remark 5
     constant 7/3. These look approachable by the Euler-integral / Hadamard-product representation (EHMM III §6) and a
     quadratic transformation linking 2F1(¼,¾;1;·) to K. This is the same machinery that produced RWZ (30)–(31) and
     our (1/2,1/2) proof.
  2. **(1/5,2/5)**. The deep input, the period identification, is proved (BGK 2025). What remains is bookkeeping:
     the F_j ↔ T_µ dictionary, the cocycle constant, and Manin's normalisation.
  3. **(1/3,1/4)**. The L-values are proved (CM). One needs Γ-evaluations of `F_k(1)` for 4F3(⅓,⅔,¼,¾;1^3;1), which
     is plausible via the CM structure of the factors 2F1(⅓,⅔) and 2F1(¼,¾).

## 4. Reproducibility (scratch, not committed)
All files below are in `/tmp/claude-0/-home-user-weighted-prime-races/91ffb489-891a-5123-b3c1-4a4206e0738c/scratchpad/l13_work/`:
* `fj2.py`: F_j(1) by exact ε-series recursion plus Richardson extrapolation. It reproduces S_0 = F_2/(2π^2) for
  (1/2,1/2) to 28 digits.
* `lv.gp`, `lv.out`: PARI L-values.
* `allF.py`, `allF.out`: ratio identification.
* `tw2.gp`, `tw_res_*.out`: twisted L-values.

Downloaded texts are in `scratchpad/l13_*/p.txt` (untrusted; read with grep/sed only).
